# -*- coding: utf-8 -*-
"""Classes and functions for manipulating/processing Tree."""

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.3'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json


class Node:
    """A node is an element of a tree structure. It is linked with its parent and its children and it is described by several metadata."""

    def __init__(self, name=None, parent_node=None, children_nodes=None, metadata=None):
        """
        Build and return an instance of Node.

        :param name: The node name.
        :type name: str
        :param parent_node: The parent node.
        :type parent_node: Node
        :param children_nodes: List of children nodes.
        :type children_nodes: list
        :param metadata: The metadata.
        :type metadata: dict
        :return: The new instance.
        :rtype: node.Node
        """
        self.children = list()
        if children_nodes is not None:
            for current_child in children_nodes:
                self.addChild(current_child)
        self.name = name
        self.metadata = metadata if metadata is not None else dict()
        self.parent = None
        if parent_node is not None:
            parent_node.addChild(self)

    def __str__(self):
        """
        Return the “informal” or nicely printable string representation of the object.

        :return: The printable string representation of the object.
        :rtype: str
        """
        node_str = "-" if self.name is None else self.name
        node_str += "\n\tParent=" + ("-" if self.parent is None or self.parent.name is None else self.parent.name)
        node_str += "\n\tChilds=" + ", ".join(self.children.keys())
        metadata = list()
        for key in sorted(self.metadata.keys()):
            metadata.append(key + ":" + str(self.metadata[key]))
        node_str += "\n\tMetadata=" + ", ".join(metadata)
        return node_str

    def hasChild(self, name=None):
        """
        Return True if the node has the specified child (name is set) or at least one child (name is None).

        :param name: The name of the searched child.
        :type name: str
        :return: True if the node has the specified child (name is set) or at least one child (name is None).
        :rtype: bool
        """
        if name is None:
            return len(self.children) > 0
        else:
            children_names = [child.name for child in self.children]
            return name in children_names

    def getChildByName(self, name):
        """
        Return the specified child node.

        :param name: Name of the searched child.
        :type name: str
        :return: The child.
        :rtype: node.Node
        """
        selected_child = None
        for child in self.children:
            if child.name is not None and child.name == name:
                selected_child = child
        if selected_child is None:
            raise Exception(
                "The node doesn't have a child named '{}'.".format(name)
            )
        return selected_child

    def getAncestors(self):
        """
        Return the ancestors of the node.

        :return: List ancestors nodes. The nodes are in parent to child order.
        :rtype: list
        """
        ancestors = list()
        if self.parent is not None:
            ancestors.extend(self.parent.getAncestors())
            ancestors.extend([self.parent])
        return ancestors

    def getDescendants(self, depth=1):
        """
        Return the node descendants with the provided depth from the node. Example: depth=1 returns all the children of the node ; depth=2 returns all the grandchildren of the node.

        :param depth: The selected depth.
        :type depth: int
        :return: Descendants.
        :rtype: list
        """
        descendants = list()
        if depth == 1:
            descendants = self.children
        elif depth > 1:
            for child in self.self.children:
                descendants.extend(child.getDescendants(depth - 1))
        return descendants

    def getLeaves(self):
        """
        Return leaves.

        :return: Leaves.
        :rtype: list
        """
        leaves = list()
        if not self.hasChild():
            leaves = [self]
        else:
            for child in self.children:
                leaves.extend(child.getLeaves())
        return leaves

    def getDepth(self):
        """
        Return the depth of the node (= the branch length / = the number of ancestors nodes).

        :return: The depth of the node. The depth for the root element is 0.
        :rtype: int
        """
        if self.parent is None:
            return 0
        else:
            return(self.parent.getDepth() + 1)

    def addChild(self, child):
        """
        Add a node as child.

        :param child: The added node.
        :type child: node.Node
        """
        if child.name is not None:  # Check duplication for named nodes
            children_names = [child.name for child in self.children]
            if child.name in children_names:
                raise Exception(
                    "Duplicated child name '{}' in node '{}'.".format(
                        child.name, self.name
                    )
                )
        child.parent = self
        self.children.append(child)

    @staticmethod
    def fromDict(tree):
        """
        Return nodes representing the tree.

        :param tree: [dict] Each node is dict with following form: {"name": NAME, "metadata":{...}, "children": OTHER_NODE}. Each attribute is optional.
        :return: The node representing the root.
        :rtype: node.Node
        """
        # Current node
        curr_node = Node()
        if "name" in tree:
            curr_node.name = tree["name"]
        if "metadata" in tree:
            for met_key, met_val in tree["metadata"].items():
                curr_node.metadata[met_key] = met_val
        # Chilren
        if "children" in tree and len(tree["children"]) > 0:
            for child in tree["children"]:
                curr_node.addChild(
                    Node.fromDict(child)
                )
        return curr_node

    @staticmethod
    def fromClusterNode(tree, id_to_name=None, distance_tag="dist", _parent_dist_from_root=0, _root_dist_from_leaves=None):
        """
        Return nodes representing the tree.

        :param tree: The tree returned by numpy.
        :type tree: scipy.cluster.hierarchy.ClusterNode
        :param id_to_name: The link between numpy's node id and the node name.
        :type id_to_name: dict
        :param distance_tag: Tag used to store distance in metadata.
        :type distance_tag: str
        :param _parent_dist_from_root: Distance between root and parent node. By default the first node provided is view as root.
        :type _parent_dist_from_root: float
        :param _root_dist_from_leaves: Maximum distance between root and leaves. By default the first node provided is view as root.
        :type _root_dist_from_leaves: float
        :return: The node representing the root.
        :rtype: node.Node
        """
        if _root_dist_from_leaves is None:  # First node is root
            _root_dist_from_leaves = tree.dist
        # Current node
        node_name = None
        if id_to_name is not None and tree.get_id() in id_to_name:
            node_name = id_to_name[tree.get_id()]
        curr_dist_from_root = _root_dist_from_leaves - tree.dist
        curr_node = Node(
            node_name,
            metadata={
                distance_tag: curr_dist_from_root - _parent_dist_from_root,
                "dist_from_root": curr_dist_from_root
            }
        )
        # Chilren
        if tree.get_left():
            curr_node.addChild(
                Node.fromClusterNode(tree.get_left(), id_to_name, distance_tag, curr_dist_from_root, _root_dist_from_leaves)
            )
        if tree.get_right():
            curr_node.addChild(
                Node.fromClusterNode(tree.get_right(), id_to_name, distance_tag, curr_dist_from_root, _root_dist_from_leaves)
            )
        return curr_node

    def toDict(self):
        """
        Return a dictionary representing the instance.

        :return: The dictionary representing the instance with following form: {"name": NAME, "metadata":{...}, "children": OTHER_NODE}.
        :rtype: dict
        """
        return {
            "name": self.name,
            "metadata": self.metadata,
            "children": [child.toDict() for child in self.children]
        }

    def toNewick(self, distance_tag=None):
        """
        Return the representation of the tree rooted by the node in newick format.

        :param distance_tag: The metadata tag for the node distance (default: 'dist'). The distance is not necessary to use this method.
        :type distance_tag: str
        :returns: Newick representation of the tree.
        :rtype: str
        """
        if distance_tag is None:
            distance_tag = "dist"
        node_newick = ""
        if self.hasChild():
            children_newick = list()
            for child in self.children:
                children_newick.append(child.toNewick())
            node_newick = '({})'.format(','.join(children_newick))
        if self.name is not None:
            node_newick += '"' + self.name + '"'
        if distance_tag in self.metadata:
            node_newick += ':{}'.format(self.metadata[distance_tag])
        return node_newick

    def toExtendedNewick(self):
        """
        Return the representation of the tree rooted by the node in extended newick format. In extended newick the distance tag is replaced by a json represetation of the metadata.

        :return: Extended newick representation of the tree.
        :rtype: str
        """
        node_newick = ""
        if self.hasChild():
            children_newick = list()
            for child in self.children:
                children_newick.append(child.toExtendedNewick())
            node_newick = '({})'.format(','.join(children_newick))
        if self.name is not None:
            node_newick += '"' + self.name + '"'
        if len(self.metadata) != 0:
            node_newick += ':{}'.format(json.dumps(self.metadata))
        return node_newick
