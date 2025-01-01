from toolbox.training import pulseDrive

tr = pulseDrive.from_config("7_1/","minimized.txt")
for _ in range(5):
    tr.pick_modified_cell()

tr.run(iterations = 10)