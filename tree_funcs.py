import tree_node


# createTree() builds the FP-tree from the data_set by filtering the genes that satisfy the
# min_sup and length conditions.
def create_tree(data_set, length, min_sup=1):
    header_table = {}  # {gene: no. of genomes it appears in}

    # this pass counts frequency of occurrence
    for genome in data_set:  # genome id
        local_set = set()
        for window in data_set[genome]:
            for gene in window:
                local_set.add(gene)
        for gene in local_set:
            header_table[gene] = header_table.get(gene, 0) + 1

    # now, header_table is {gene: no. of genomes it appears in}

    for gene in list(header_table):  # remove items not meeting minSup
        if header_table[gene] < min_sup:
            del (header_table[gene])
    freq_gene_set = set(header_table.keys())

    if len(freq_gene_set) == 0:
        return None  # if no genes meet min support --> get out

    ret_tree = tree_node.TreeNode('Null Set', 1, None)  # create tree

    for genome, windows in data_set.items():
        for window in windows:
            filtered_window_d = {}  # {gene: no of genomes it's in}
            for gene in window:
                if gene in freq_gene_set:
                    filtered_window_d[gene] = header_table[gene]
                    # this makes sure that each gene appears only once in ordered items

            if len(filtered_window_d) > length:  # there are frequent genes in this window - at least l
                ordered_path = [v[0] for v in sorted(filtered_window_d.items(), key=lambda p: p[1], reverse=True)]
                # (ordered_path is the window that is ordered in a descending order based on
                # frequency of each gene in the genomes)
                update_tree(ordered_path, ret_tree, genome)  # populate tree with ordered_path

    return ret_tree


# updateTree() grow the Fp-tree with a window of genes.
def update_tree(window, in_tree, genome_num):
    if window[0] in in_tree.children:  # check if window[0] in retTree.children
        if in_tree.children[window[0]].last_update != genome_num:
            in_tree.children[window[0]].inc(1)  # increase count by one
            in_tree.children[window[0]].last_update = genome_num
    else:  # add window[0] to in_tree.children
        in_tree.children[window[0]] = tree_node.TreeNode(window[0], 1, in_tree)
        in_tree.children[window[0]].last_update = genome_num

    if len(window) > 1:  # call updateTree() with remaining ordered items
        update_tree(window[1::], in_tree.children[window[0]], genome_num)


# traverse the tree with the min_sup condition:
def dfs(node, final_ans, curr_path, min_sup):
    if node.parent is not None:
        curr_path.append(node.name)

    if len(node.children) == 0:  # leaf
        if node.count >= min_sup:  # here we apply the min_sup constraint to a path(hitchhiker's group)
            final_ans.append([curr_path, node.count])

        return

    for child in node.children.values():
        new_path = curr_path[:]
        dfs(child, final_ans, new_path, min_sup)
