def printl(list_, sep=' '):
    print(sep.join([str(i) for i in list_]))


def print_list_as_path(euler_cycle):
    print('->'.join(euler_cycle))