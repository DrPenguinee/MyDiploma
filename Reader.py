with open('vneut1_sample.dat', 'r') as f:
    i = 0
    for i in range(100):
        if i > 100:
            break
        item = f.readline().split()
        print(item)
        i += 1