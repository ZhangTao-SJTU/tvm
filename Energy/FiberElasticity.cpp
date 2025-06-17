/* ---------------------------------------------------------------------------------
 * Copyright 2021-2023 Tao Zhang
 *
 * This file is part of TVM.
 *
 * TVM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation,
 * either version 3 of the License,
 * or (at your option) any later version.
 *
 * TVM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TVM. If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
 * Coauthor: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu
 * ---------------------------------------------------------------------------------
 */

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_map>

#include "FiberElasticity.h"

using namespace std;

FiberElasticity::FiberElasticity(Run * run) {
    run_ = run;
    ks_ = 0.;
    l0_ = 0.;
    kb_ = 0.;
}

int     FiberElasticity::updateForces() {
    // reset all stretchForce values in nodes
    for (auto node : run_->nodes_) {
        for (int m = 0; m < 3; m++) {
            node->stretchForce_[m] = 0.;
        }
    }
    // reset all bendForce values in nodes
    for (auto node : run_->nodes_) {
        for (int m = 0; m < 3; m++) {
            node->bendForce_[m] = 0.;
        }
    }

    // update stretchForce values
    for (auto fiber : run_->fibers_) {
        for (int i = 0; i < fiber->nodes_.size() - 1; i++) {
            Node * ni = fiber->nodes_[i];
            Node * nj = fiber->nodes_[i + 1];
            updateStretchForces(ni, nj);
        }
    }

    // update bendForce values
    for (auto fiber : run_->fibers_) {
        for (int i = 1; i < fiber->nodes_.size() - 1; i++) {
            Node * nj = fiber->nodes_[i - 1];
            Node * ni = fiber->nodes_[i];
            Node * nk = fiber->nodes_[i + 1];
            updateBendForces(nj, ni, nk);
        }
    }

    return 0;
}

int FiberElasticity::updateStretchForces(Node * n0, Node * n1) {
    double dx[3];
    dx[0] = n1->position_[0] - n0->position_[0];
    dx[1] = n1->position_[1] - n0->position_[1];
    dx[2] = n1->position_[2] - n0->position_[2];
    run_->box_->resetDistance(dx);
    double l = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    double Fabs = ks_ * (-1. + l0_ / l);
    n0->stretchForce_[0] = n0->stretchForce_[0] - Fabs * dx[0];
    n1->stretchForce_[0] = n1->stretchForce_[0] + Fabs * dx[0];
    n0->stretchForce_[1] = n0->stretchForce_[1] - Fabs * dx[1];
    n1->stretchForce_[1] = n1->stretchForce_[1] + Fabs * dx[1];
    n0->stretchForce_[2] = n0->stretchForce_[2] - Fabs * dx[2];
    n1->stretchForce_[2] = n1->stretchForce_[2] + Fabs * dx[2];

    return 0;
}

int FiberElasticity::updateBendForces(Node * nj, Node * ni, Node * nk) {
    double dxD[3];
    double dxU[3];
    for (int m = 0; m < 3; m++) {
        dxD[m] = nj->position_[m] - ni->position_[m];
        dxU[m] = nk->position_[m] - ni->position_[m];
    }
    run_->box_->resetDistance(dxD);
    run_->box_->resetDistance(dxU);

    double lD = sqrt(dxD[0] * dxD[0] + dxD[1] * dxD[1] + dxD[2] * dxD[2]);
    double lU = sqrt(dxU[0] * dxU[0] + dxU[1] * dxU[1] + dxU[2] * dxU[2]);
    double dP = dxD[0] * dxU[0] + dxD[1] * dxU[1] + dxD[2] * dxU[2];
    double cosT = dP / (lD * lU);

    double FbendC = -1. * kb_ * (cosT + (double)1.);

    double fbUX = (cosT * dxU[0] / lU - dxD[0] / lD) * FbendC / lU;
    double fbUY = (cosT * dxU[1] / lU - dxD[1] / lD) * FbendC / lU;
    double fbUZ = (cosT * dxU[2] / lU - dxD[2] / lD) * FbendC / lU;
    double fbDX = (cosT * dxD[0] / lD - dxU[0] / lU) * FbendC / lD;
    double fbDY = (cosT * dxD[1] / lD - dxU[1] / lU) * FbendC / lD;
    double fbDZ = (cosT * dxD[2] / lD - dxU[2] / lU) * FbendC / lD;

    ni->bendForce_[0] = ni->bendForce_[0] + (fbUX + fbDX);
    ni->bendForce_[1] = ni->bendForce_[1] + (fbUY + fbDY);
    ni->bendForce_[2] = ni->bendForce_[2] + (fbUZ + fbDZ);

    nk->bendForce_[0] = nk->bendForce_[0] - fbUX;
    nk->bendForce_[1] = nk->bendForce_[1] - fbUY;
    nk->bendForce_[2] = nk->bendForce_[2] - fbUZ;

    nj->bendForce_[0] = nj->bendForce_[0] - fbDX;
    nj->bendForce_[1] = nj->bendForce_[1] - fbDY;
    nj->bendForce_[2] = nj->bendForce_[2] - fbDZ;

    return 0;
}

int FiberElasticity::updateEnergy() {
    stretchEnergy_ = 0.;

    for (auto fiber : run_->fibers_) {
        for (int i = 0; i < fiber->nodes_.size() - 1; i++) {
            Node * ni = fiber->nodes_[i];
            Node * nj = fiber->nodes_[i + 1];
            double dx[3];
            dx[0] = nj->position_[0] - ni->position_[0];
            dx[1] = nj->position_[1] - ni->position_[1];
            dx[2] = nj->position_[2] - ni->position_[2];
            run_->box_->resetDistance(dx);
            double l = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
            stretchEnergy_ += 0.5 * ks_ * pow(l - l0_, 2.0);
        }
    }

    bendEnergy_ = 0.;
    for (auto fiber : run_->fibers_) {
        for (int i = 1; i < fiber->nodes_.size() - 1; i++) {
            Node * nj = fiber->nodes_[i - 1];
            Node * ni = fiber->nodes_[i];
            Node * nk = fiber->nodes_[i + 1];
            double dxD[3];
            double dxU[3];
            for (int m = 0; m < 3; m++) {
                dxD[m] = nj->position_[m] - ni->position_[m];
                dxU[m] = nk->position_[m] - ni->position_[m];
            }
            run_->box_->resetDistance(dxD);
            run_->box_->resetDistance(dxU);
            double lD = sqrt(dxD[0] * dxD[0] + dxD[1] * dxD[1] + dxD[2] * dxD[2]);
            double lU = sqrt(dxU[0] * dxU[0] + dxU[1] * dxU[1] + dxU[2] * dxU[2]);
            double dP = dxD[0] * dxU[0] + dxD[1] * dxU[1] + dxD[2] * dxU[2];
            double cosT = dP / (lD * lU);
            bendEnergy_ += 0.5 * kb_ * pow(cosT + 1.0, 2.0);
        }
    }

    return 0;
}