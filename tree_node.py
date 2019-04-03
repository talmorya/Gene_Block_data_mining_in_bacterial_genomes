class TreeNode:

    def __init__(self, name_value, num_occur, parent_node):
        self.name = name_value  # gene cog number
        self.count = num_occur  # in how many genomes this gene appeared in a valid window
        self.parent = parent_node
        self.children = {}
        self.last_update = 0  # the last genome that updated this node in the creation of the tree

    # increments the count variable with a given amount
    def inc(self, num_occur):
        self.count += num_occur

    # display tree in text. Useful for debugging.
    def disp(self, ind=1):
        print('  ' * ind, self.name, ' ', self.count)
        for child in self.children.values():
            child.disp(ind + 1)
