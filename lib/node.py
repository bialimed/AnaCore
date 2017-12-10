#
# Copyright (C) 2015 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frederic.escudie@iuct-oncopole.fr'
__status__ = 'prod'


import json


class Node:
    """
    @summary: A node is an element of a tree structure. It is linked with its
    parent and its children and it is described by several metadata.
    """
    def __init__(self, name, parent_node=None, children_nodes=None, metadata=None):
        """
        @param name: [str] The node name.
        @param parent_node: [Node] The parent node.
        @param children_nodes: [list] List of children nodes.
        @param metadata: [dict] The metadata.
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
        @summary: Returns a string representation of the node.
        @return: [str] The representation of the node.
        """
        node_str = self.name
        node_str += "\n\tParent=" + (self.parent.name if self.parent is not None else "-")
        node_str += "\n\tChilds=" + ", ".join(self.children.keys())
        metadata = list()
        for key in sorted(self.metadata.keys()):
            metadata.append(key + ":" + str(self.metadata[key]))
        node_str += "\n\tMetadata=" + ", ".join(metadata)
        return node_str

    def hasChild(self, name=None):
        """
        @summary: Returns true if the node has the specified child (name is
        set) or at least one child (name is None).
        @param name: [str] the name of the searched child.
        @return: [bool]
        """
        if name is None:
            return len(self.children) > 0
        else:
            children_names = [child.name for child in self.children]
            return name in children_names

    def getChildByName(self, name):
        """
        @summary: Returns the specified child node.
        @param name: [str] the name of the searched child.
        @return: [Node] The child.
        """
        selected_child = None
        for child in self.children:
            if child.name is not none and child.name == name:
                selected_child = child
        if selected_child is None:
            raise Exception(
                "The node doesn't have a child named '{}'.".format(name)
            )
        return selected_child

    def getAncestors(self):
        """
        @summary: Returns the ancestors of the node.
        @return: [list] The list ancestors nodes. The nodes are in parent to
        child order.
        """
        ancestors = list()
        if self.parent is not None:
            ancestors.extend(self.parent.getAncestors())
            ancestors.extend([self.parent])
        return ancestors

    def getDescendants(self, depth=1):
        """
        @summary: Returns the node descendants with the provided depth from the
        node. Example: depth=1 returns all the children of the node ; depth=2
        returns all the grandchildren of the node.
        @param: [int] The selected depth.
        @return: [list] The nodes of descendants.
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
        @summary: Returns leaves.
        @return: [list] The nodes of leaves.
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
        @summary: Returns the depth of the node (= the branch length/ = the
        number of ancestors nodes).
        @return: [int] The depth of the node. The depth for the root element is 0.
        """
        if self.parent is None:
            return 0
        else:
            return(self.parent.getDepth() + 1)

    def addChild(self, child):
        """
        @summary: Adds a node as child.
        @param child: [Node] The added node.
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
        @summary: Returns nodes representing the tree.
        @param tree: [dict] Each node is dict with following form: {"name":
        NAME, "metadata":{...}, "children": OTHER_NODE}. Each attribute is
        optional.
        @return: [Node] The node representing the root.
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
            for child in children:
                curr_node.addchild(
                    Node.fromDict(child)
                )
        return curr_node

    @staticmethod
    def fromClusterNode(tree, id_to_name=None, distance_tag="dist", _metadata=None):
        """
        @summary: Returns nodes representing the tree.
        @param tree: [ClusterNode] The tree returned by numpy.
        @param id_to_name: [dict] The link between numpy's node id and the node
        name.
        @param distance_tag: [str] Tag used to store distance in metadata.
        @param _metadata: [dict] Intern param used to pass information between
        node and children.
        @return: [Node] The node representing the root.
        """
        # Current node
        node_name = None
        if id_to_name is not None and tree.get_id() in id_to_name:
            node_name = id_to_name[tree.get_id()]
        curr_node = Node(node_name, metadata=_metadata)
        # Chilren
        if tree.get_left():
            curr_node.addChild(
                Node.fromClusterNode(tree.get_left(), id_to_name, distance_tag, {distance_tag: tree.dist})
            )
        if tree.get_right():
            curr_node.addChild(
                Node.fromClusterNode(tree.get_right(), id_to_name, distance_tag, {distance_tag: tree.dist})
            )
        return curr_node

    def toDict(self):
        """
        @summary: Returns a dictionary representing the instance.
        @return: [dict] The dictionary representing the instance with following
        form: {"name": NAME, "metadata":{...}, "children": OTHER_NODE}.
        """
        return {
            "name": self.name,
            "metadata": self.metadata,
            "children": [child.toDict() for child in self.children]
        }

    def toNewick(self, distance_tag=None):
        """
        @summary: Returns the representation of the tree rooted by the node in
        newick format.
        @param distance_tag: [str] The metadata tag for the node distance
        (default: 'dist'). The distance is not necessary to use this method.
        @returns: [str] the newick representation of the tree.
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
        @summary: Returns the representation of the tree rooted by the node in
        extended newick format. In extended newick the distance tag is replaced
        by a json represetation of the metadata.
        @returns: [str] the extended newick representation of the tree.
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
