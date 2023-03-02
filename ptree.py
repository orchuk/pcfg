class Node:
    def __init__(self, key, children):
        self.key = key
        self.children = children

    def __repr__(self):
        if self.children:
            return '[{} {}]'.format(self.key, ' '.join(map(repr, self.children)))
        else:
            return self.key


class PTree:
    def __init__(self, root=None, probability=0):
        self.root = root
        self.probability = probability

    def __repr__(self):
        return '({}): {}'.format(self.probability, repr(self.root))
