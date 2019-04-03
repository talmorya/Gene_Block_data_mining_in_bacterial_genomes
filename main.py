import sys

import cog_funcs
import pre_process
import tree_funcs
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # get the arguments from the user:
    d = int(sys.argv[1])
    length = int(sys.argv[2])  # l
    min_sup = int(sys.argv[3])
    query = sys.argv[4:]

    # get the data_set to work with from the pre-processing:
    data_set = {}  # {window: genome id}
    data_set = pre_process.start(d, query, length)

    # create the tree from the data_set:
    tree = tree_funcs.create_tree(data_set, length - len(query), min_sup)
    if tree is None:  # if the tree is empty:
        print("there are no frequent hitchhiker's")
        sys.exit()

    # get the hitchhikers paths from the tree:
    freq_paths = []
    tree_funcs.dfs(tree, freq_paths, [], min_sup)
    if len(freq_paths) == 0:  # if there isn't any that appears at least min_sup
        print("there are no frequent hitchhiker's")
        sys.exit()

    # top_paths are the paths ordered by the frequency of the paths
    top_paths = sorted(freq_paths, key=lambda p: p[1], reverse=True)

    # print the info of the top three paths:
    top_three = set(query)
    for path in top_paths[:3]:
        for gene in path[0]:
            top_three.add(gene)

    # a dictionary of the info of the cogs in the top three paths:
    cogs_info = cog_funcs.get_cogs_info(top_three)

    print("The query is: ")
    for q in query:
        print(q + ": " + '. '.join(cogs_info[q]))

    i = 1
    for path in top_paths[:3]:
        print("*********************************************************")
        print("Result no. ", i, " : ", path[0])
        print("The number of genomes containing this sequence is: ", path[1], '\n')
        for gene in path[0]:
            print(gene + ": " + '. '.join(cogs_info[gene]))

        i += 1

    # plot the graph:
    x = []
    y = []
    for path in freq_paths:
        x.append('\n'.join(path[0]))
        y.append(int(path[1]))

    plt.xlabel('Hitchhikers sequence')
    plt.ylabel('no. of genomes containing')

    plt.bar(x, y)
    for a, b in zip(x, y):
        plt.text(a, b, str(b))
    plt.show()
