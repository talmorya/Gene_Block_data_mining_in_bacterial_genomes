data_set = {}
d = 0
l = 0
query = []


def start(window_size, q, length):
    global d, query, l
    d = window_size
    l = length
    query = q

    # open files and process them:
    open_files('/Users/Shani/PycharmProjects/Gene_Block_data_mining_in_bacterial_genomes/cog_words_bac.txt')
    # open_files('/Users/Shani/PycharmProjects/Gene_Block_data_mining_in_bacterial_genomes/cog_words_plasmid.txt')

    return data_set


def open_files(path):

    with open(path) as f:
        # use readline() to read the first line
        line = f.readline()

        # while there is a next line parse it
        while line:
            parse_line(line)
            line = f.readline()


def parse_line(line):
    line_elements = line.split('\t')
    genome_num = line_elements[0].split('#')[-1]
    if line_elements[-1] == '\n':
        word = line_elements[1:-1]
    else:
        word = line_elements[1:]

    parse_word(word, genome_num)


def parse_word(word, genome_num):
    # if len of word is smaller than the query it's irrelevant
    if len(word) < len(query):
        return

    # global d , query
    q_dic = {}
    for g in query:
        q_dic[g] = 0

    first = word[0]
    last_window = False

    # go over first window
    for i in range(d):

        if i >= len(word)-1:  # if word< 1st window
            if 0 in q_dic.values():  # non valid window, do nothing
                return
            else:
                add_window(word, genome_num)  # add the window if it is valid
                return

        if word[i] in q_dic:
            q_dic[word[i]] += 1

    next_l = word[d]

    # if the first window is valid insert it to dataSet
    if 0 not in q_dic.values():
        add_window(word[:d], genome_num)
        last_window = True

    idx = d
    while idx < len(word):
        if last_window:
            if first in q_dic.keys():
                q_dic[first] -= 1
                if q_dic[first] == 0:
                    if first == next_l:
                        add_window(word[idx - d + 1:idx + 1], genome_num)
                        last_window = True
                    else:
                        last_window = False
                else:
                    add_window(word[idx - d + 1:idx + 1], genome_num)
                    last_window = True
            else:  # first of last window is irrelevant, good window
                add_window(word[idx - d + 1:idx + 1], genome_num)
                last_window = True

        else:  # last window was false:
            if list(q_dic.values()).count(0) == 1 and next_l in q_dic and q_dic[next_l] == 0:
                add_window(word[idx - d + 1:idx + 1], genome_num)
                last_window = True
            else:
                last_window = False

        # updating shit for next window
        idx += 1
        if idx < len(word):  # not last round
            if next_l in q_dic.keys():
                q_dic[next_l] += 1
            first = word[idx - d]
            next_l = word[idx]


def add_window(window, genome_num):
    global data_set
    final_window = list(filter(lambda x: (x != 'X') and (x not in query), window))

    # we don't need to add the paths that are smaller than l to the tree
    if len(final_window) < l-len(query):
        return

    final_window.sort()

    if genome_num in data_set.keys():
        if final_window not in data_set[genome_num]:
            data_set[genome_num].append(final_window)
    else:
        data_set[genome_num] = [final_window]
