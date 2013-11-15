class Node:
    """
    Tree node: left and right child + data which can be any object
    """
    def __init__(self, data):
        """
        Node constructor, data is a two element tuple

        @param data node data object
        """
        self.left = None
        self.right = None
        self.data = data
        self.value = max()

    # def inRegion(self,seg1,seg2):
    #      ''' evaluation whether seg1 in in seg2: seg=chr_start_end'''
    #      seg1Array = map(int,seg1.split("_"))
    #      seg2Array = map(int,seg2.split("_"))
    #      if len(seg1Array) == len(seg2Array) and len(seg1Array) == 3:
    #        if seg1Array[0] == seg2Array[0]:
    #             if seg1Array[2] < seg2Array[1] or seg1Array[1] > seg2Array[2]:
    #                 flag = -1
    #             elif seg1Array[1] == seg2Array[1] and seg1Array[2]  == seg2Array[2]:
    #                 flag = 0 
    #             else:
    #                 flag = 1
    #        else:
    #           flag = 0
    #      else:
    #        flag = 0
    #      return(flag)


    def insert(self, data):
        """
        Insert new node with data

        @param data node data object to insert
        """
        if data[0] :
            if self.left is None:
                self.left = Node(data)
            else:
                self.left.insert(data)
        else:
            if self.right is None:
                self.right = Node(data)
            else:
                self.right.insert(data)

    def lookup(self, data, parent=None):
        """
        Lookup node containing data

        @param data node data object to look up
        @param parent node's parent
        @returns node and node's parent if found or None, None
        """
        if self.inRegion(data,self.data) < 0:
            if self.left is None:
                return None, None
            return self.left.lookup(data, self)
        elif self.inRegion(data,self.data) > 0:
            if self.right is None:
                return None, None
            return self.right.lookup(data, self)
        else:
            return self, parent

    def delete(self, data):
        """
        Delete node containing data

        @param data node's content to delete
        """
        # get node containing data
        node, parent = self.lookup(data)
        if node is not None:
            children_count = node.children_count()
            if children_count == 0:
                # if node has no children, just remove it
                # check if it is not the root node
                if parent.left is node:
                    parent.left = None
                else:
                    parent.right = None
                del node
            elif children_count == 1:
                # if node has 1 child
                # replace node by its child
                if node.left:
                    n = node.left
                else:
                    n = node.right
                if parent.left is node:
                    parent.left = n
                else:
                    parent.right = n
                del node
            else:
                # if node has 2 children
                # find its successor
                parent = node
                successor = node.right
                while successor.left:
                    parent = successor
                    successor = successor.left
                # replace node data by its successor data
                node.data = successor.data
                # fix successor's parent node child
                if parent.left == successor:
                    parent.left = successor.right
                else:
                    parent.right = successor.right

    def compare_trees(self, node):
        """
        Compare 2 trees

        @param node tree to compare
        @returns True if the tree passed is identical to this tree
        """
        if node is None:
            return False
        if self.data != node.data:
            return False
        res = True
        if self.left is None:
            if node.left:
                return False
        else:
            res = self.left.compare_trees(node.left)
        if self.right is None:
            if node.right:
                return False
        else:
            res = self.right.compare_trees(node.right)
        return res
                
    def print_tree(self):
        """
        Print tree content inorder
        """
        if self.left:
            self.left.print_tree()
        print self.data,
        if self.right:
            self.right.print_tree()

    def tree_data(self):
        """
        Generator to get the tree nodes data
        """
        # we use a stack to traverse the tree in a non-recursive way
        stack = []
        node = self
        while stack or node: 
            if node:
                stack.append(node)
                node = node.left
            else: # we are returning so we pop the node and we yield it
                node = stack.pop()
                yield node.data
                node = node.right

    def children_count(self):
        """
        Return the number of children

        @returns number of children: 0, 1, 2
        """
        cnt = 0
        if self.left:
            cnt += 1
        if self.right:
            cnt += 1
        return cnt
