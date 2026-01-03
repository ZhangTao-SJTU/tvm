from scripts.diagnostic_plots import plot_all


def main():
    run_dir = "tests/test_multiple_patterns_4_cells_v2/"
    n_cells = 4
    out_dir = run_dir + "plots/"
    plot_all(run_dir, out_dir, n_cells)
if __name__ == "__main__":
    main()