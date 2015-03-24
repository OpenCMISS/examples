""" Module for reading exnode and exelem files used by cmgui

    Has only been used for reading exfiles produced by OpenCMISS
    and doesn't attempt to be able to read any valid file.
"""

import gzip
import numpy as np
import re

__all__ = [
        "Exregion"  
        "Exnode",
        "ExnodeField",
        "ExnodeComponent",
        "ExnodeNode",
        "Exelem",
        "ExelemField",
        "ExelemComponent",
        "ExelemElement",
        "ExfileError",
        ]

class Exregion(object):
    """ Store and retrieve data from an exelem file
    """
    def __init__(self, filepath):
        self.fields = []
        self.elements = []
        self.sections = []
        with FileWithLineNumber(filepath, 'r') as f:
            self._read_header(f)
            self.num_element_values = self._calc_num_element_values()
            while True:
                try:
                    self.sections.append(ExnodeSection(f, self))
                except ExfileError:
                    break
# Go back to the past 2 lines    
            f.rollbacktwice()
# read the exelem header data                
            self.num_dims = read_regex(f, r'Shape.\s+Dimension=\s*([0-9]+)')
            self.num_scale_factor_sets = int(read_regex(f, r'#Scale factor sets=\s*([0-9]+)'))
            self.num_scale_factors = 0
            for i in range(self.num_scale_factor_sets):
                self.num_scale_factors += int(read_regex(f, r'#Scale factors\s*=\s*([0-9]+)'))
            self.num_nodes = int(read_regex(f, r'#Nodes=\s*([0-9]+)'))
            self.num_fields = int(read_regex(f, r'#Fields=\s*([0-9]+)'))
            for i in range(self.num_fields):
                field = ExelemField(f, self)
                self.fields.append(field)
            while True:
                try:
                    self._read_element(f)
                except EOFError:
                    break
        self.num_elements = len(self.elements)
        self.num_nodes = sum(section.num_nodes for section in self.sections)

    def _read_header(self, f):
        line = f.readline().strip()
        regex = r'Group name|Region: ([A-Za-z0-9_\:\.\/]+)$'
        match = re.search(regex, line)
#Skip root region specification        
        while match is None:
            line = f.readline().strip()
            match = re.search(regex, line)
        if len(match.groups()) == 1:
            self.group_name = match.group(1)
        else:
            self.group_name = match.groups()
#ZINC generated exregion files have  !#nodeset nodes, check for that
        line = f.readline().strip()
        regex = r'nodeset nodes$'
        match = re.search(regex, line)
        if match is None:
            f.rollback();            
            
    def node_values(self, field_name, component_name, node_num):
        """ Return all the field component derivative values
            at the given node number
        """
        for section in self.sections:
            try:
                return section.node_values(field_name, component_name, node_num)
            except NodeNotFound:
                pass
        raise ValueError("Node %d not found in any exnode section." % node_num)

    def node_value(self, field_name, component_name, node_num,
            derivative_number=1):
        """ Return the field component value at the given node and derivative
            Derivatives are numbered from 1, with 1 being no derivative.
        """
        for section in self.sections:
            try:
                return section.node_value(field_name, component_name, node_num,
                        derivative_number)
            except NodeNotFound:
                pass
        raise ValueError("Node %d not found in any exnode section." % node_num)
            

    def _calc_num_element_values(self):
        num_values = 0
        for field in self.fields:
            for component in field.components:
                if component.component_type == 'grid based':
                    component_num_values = np.product(
                            [i + 1 for i in component.divisions])
                    num_values += component_num_values
        return num_values

    def _read_element(self, f):
        element_line = f.readline()
        if element_line == "":
            raise EOFError
        indices = map(int, element_line.split(':')[1].split())
        if indices[1] == 0 or indices[2] == 0:
            #raise ExfileError(f, "Face or line elements not supported")
            values = []
            if self.num_element_values > 0:
                expect_line(f, "Values:")
                while len(values) < self.num_element_values:
                    line = f.readline()
                    values.extend(map(float, line.split()))
    
#Ignore faces
            nodes = []   
            element_line = f.readline()
            regex = r'Faces:'
            if re.search(regex, element_line) is not None: 
                regex = r'Nodes:'
                element_line = f.readline()
                match = re.search(regex, element_line)
                while match is None:
                    element_line = f.readline()
                    match = re.search(regex, element_line)

                while len(nodes) < self.num_nodes:
                    line = f.readline()
                    nodes.extend(map(int, line.split()))
            else:
                regex = r'Nodes:'
                if re.search(regex, element_line) is not None: 
                    while len(nodes) < self.num_nodes:
                        line = f.readline()
                        nodes.extend(map(int, line.split()))
                else:
                    raise ExfileError(f, "Expected '%s', got '%s'" % ("Nodes:", line))
    
            scale_factors = []
            element_line = f.readline()
            regex = r'Scale factors:'
            if re.search(regex, element_line) is not None: 
                while len(scale_factors) < self.num_scale_factors:
                    line = f.readline()
                    scale_factors.extend(map(float, line.split()))
            else:
                f.rollback()
        
            self.elements.append(
                        ExelemElement(indices, nodes, values, scale_factors))

    def element_values(self, field_name, component_name, element_num):
        """Return the all field component values at the given element number
        """
        element = self.elements[element_num - 1]
        value_index = 0
        for field in self.fields:
            for component in field.components:
                if component.component_type == 'grid based':
                    component_num_values = np.product(
                            [i + 1 for i in component.divisions])
                    if (field.name == field_name and
                            component.name == str(component_name)):
                        return element.values[value_index:value_index + component_num_values]
                    else:
                        value_index += component_num_values
        raise ValueError("Couldn't find field and component values")




class Exnode(object):
    """ Store and retrieve data from an exnode file
    """
    def __init__(self, filepath):
        self.sections = []
        with FileWithLineNumber(filepath, 'r') as f:
            self._read_header(f)
            while True:
                try:
                    self.sections.append(ExnodeSection(f, self))
                except EOFError:
                    break
        self.num_nodes = sum(section.num_nodes for section in self.sections)

    def _read_header(self, f):
        self.group_name = read_regex(f, r'Group name|Region: ([A-Za-z0-9_\:\.]+)$')

    def node_values(self, field_name, component_name, node_num):
        """ Return all the field component derivative values
            at the given node number
        """
        for section in self.sections:
            try:
                return section.node_values(field_name, component_name, node_num)
            except NodeNotFound:
                pass
        raise ValueError("Node %d not found in any exnode section." % node_num)

    def node_value(self, field_name, component_name, node_num,
            derivative_number=1):
        """ Return the field component value at the given node and derivative
            Derivatives are numbered from 1, with 1 being no derivative.
        """
        for section in self.sections:
            try:
                return section.node_value(field_name, component_name, node_num,
                        derivative_number)
            except NodeNotFound:
                pass
        raise ValueError("Node %d not found in any exnode section." % node_num)


class ExnodeSection(object):
    def __init__(self, f, exnode):
        self.exnode = exnode
        self.fields = []
        self.nodes = []
        self._read_section_header(f)
        self.num_node_values = self._calc_num_node_values()
        while True:
            try:
                self._read_node(f)
            except EOFError:
                break
            except ExfileError:
                # Could be start of new section
                f.rollback()
                break
        self.num_nodes = len(self.nodes)

    def _read_section_header(self, f):
        try:
            self.num_fields = int(read_regex(f,
                    r'#Fields=\s*([0-9]+)'))
        except ExfileError:
            if f.readline() == '':
                raise EOFError
            else:
                raise
        for i in range(self.num_fields):
            field = ExnodeField(f, self)
            self.fields.append(field)

    def _calc_num_node_values(self):
        num_values = 0
        for field in self.fields:
            for component in field.components:
                num_values += 1 + component.num_derivatives
        return num_values

    def _read_node(self, f):
        line = f.readline().strip()
        if line == "":
            raise EOFError
        number = int(read_string_regex(f, line,
                r'Node:\s*([0-9]+)'))
        read = 0
        values = np.empty(self.num_node_values)
        while read < self.num_node_values:
            line = f.readline()
            try:
                new_values = map(float, line.split())
            except ValueError:
                raise ExfileError(f, "Expecting node values, got: %s" % line.strip())
            if read + len(new_values) > self.num_node_values:
                raise ExfileError(f, "Got more node values than expected.")
            values[read:read + len(new_values)] = new_values
            read += len(new_values)

        self.nodes.append(ExnodeNode(number, values))

    def _get_field_component(self, field_name, component_name):
        try:
            field = next(f for f in self.fields if f.name == field_name)
        except StopIteration:
            raise ValueError("Couldn't find field %s" % field_name)
        try:
            # Accept integer for component name
            component_name = str(component_name)
            component = next(
                    c for c in field.components if c.name == component_name)
        except StopIteration:
            raise ValueError("Couldn't find component %s" % component_name)
        return field, component

    def node_values(self, field_name, component_name, node_num):
        """ Return all the field component derivative values
            at the given node number
        """
        try:
            node = next(n for n in self.nodes if n.number == node_num)
        except StopIteration:
            raise NodeNotFound()
        field, component = self._get_field_component(
                field_name, component_name)

        value_index = component.value_index - 1
        component_num_values = 1 + component.num_derivatives
        return node.values[value_index:value_index + component_num_values]

    def node_value(self, field_name, component_name, node_num, derivative_number=1):
        """ Return the field component value at the given node and derivative
            Derivatives are numbered from 1, with 1 being no derivative.
        """
        try:
            node = next(n for n in self.nodes if n.number == node_num)
        except StopIteration:
            raise NodeNotFound()
        field, component = self._get_field_component(field_name, component_name)

        value_index = component.value_index - 1
        value_index += derivative_number - 1
        if not 0 < derivative_number <= component.num_derivatives + 1:
            raise ValueError("Invalid derivative number: %d" % derivative_number)
        return node.values[value_index]


class ExnodeField(object):
    """ A field definition from an exnode file
    """
    def __init__(self, f, exnode):
        self.exnode = exnode

        # Eg:
        # 1) Geometry, coordinate, rectangular cartesian, #Components=3
        self.components = []

        declaration = f.readline().split(',')
        index, self.name = read_string_regex(f, declaration[0],
                r'([0-9]+)\)\s+([A-Za-z0-9_\s\/]+)')
        self.index = int(index)
        self.num_components = int(read_string_regex(f, declaration[3],
                r'#Components\s*=\s*([0-9]+)'))
        for i in range(self.num_components):
            component = ExnodeComponent(f, self)
            self.components.append(component)

    def __repr__(self):
        return '<ExnodeField: %d. "%s", #Components=%d>' % (
                self.index, self.name, self.num_components)


class ExnodeComponent(object):
    """A field component definition from an exnode file
    """
    def __init__(self, f, field):
        self.field = field
        # Eg:
        # x.  Value index= 1, #Derivatives= 0
        # or
        # 1.  Value index= 49, #Derivatives= 7(d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)
        declaration = f.readline().strip().split(',', 1)
        self.name = read_string_regex(f, declaration[0],
                r'^([a-zA-Z0-9]+)\.')
        self.value_index = int(read_string_regex(f, declaration[0],
                r'Value index\s*=\s*([0-9]+)'))
        self.num_derivatives = int(read_string_regex(f, declaration[1],
                r'#Derivatives\s*=\s*([0-9]+)'))
        if self.num_derivatives > 0:
            derivative_name_list = read_string_regex(f, declaration[1],
                    r'#Derivatives\s*=\s*[0-9]+\s*\(([a-zA-Z0-9\,\/\s]+)\)')
            self.derivative_names = [
                    n.strip() for n in derivative_name_list.split(',')]
        else:
            self.derivative_names = []

    def __repr__(self):
        return '<ExnodeComponent: "%s". #Derivatives=%d>' % (
                self.name, self.num_derivatives)


class ExnodeNode(object):
    """A node from an exnode file
    """
    def __init__(self, number, values):
        self.number = number
        self.values = values

    def __repr__(self):
        return "<ExnodeNode %d>" % self.number


class Exelem(object):
    """ Store and retrieve data from an exelem file
    """
    def __init__(self, filepath):
        self.fields = []
        self.elements = []
        with FileWithLineNumber(filepath, 'r') as f:
            self._read_header(f)
            self.num_element_values = self._calc_num_element_values()
            while True:
                try:
                    self._read_element(f)
                except EOFError:
                    break
        self.num_elements = len(self.elements)

    def _read_header(self, f):
        self.group_name = read_regex(f,
                r'Group name|Region: ([A-Za-z0-9_\:\.]+)')
        self.num_dims = int(read_regex(f,
                r'Shape.\s+Dimension=\s*([0-9]+)'))
        self.num_scale_factor_sets = int(read_regex(f,
                r'#Scale factor sets=\s*([0-9]+)'))
        self.num_scale_factors = 0
        for i in range(self.num_scale_factor_sets):
            self.num_scale_factors += int(read_regex(f,
                r'#Scale factors\s*=\s*([0-9]+)'))
        self.num_nodes = int(read_regex(f,
                r'#Nodes=\s*([0-9]+)'))
        self.num_fields = int(read_regex(f,
                r'#Fields=\s*([0-9]+)'))
        for i in range(self.num_fields):
            field = ExelemField(f, self)
            self.fields.append(field)

    def _calc_num_element_values(self):
        num_values = 0
        for field in self.fields:
            for component in field.components:
                if component.component_type == 'grid based':
                    component_num_values = np.product(
                            [i + 1 for i in component.divisions])
                    num_values += component_num_values
        return num_values

    def _read_element(self, f):
        element_line = f.readline()
        if element_line == "":
            raise EOFError
        indices = map(int, element_line.split(':')[1].split())
        if indices[1] == 0 and indices[2] == 0:
            # raise ExfileError(f, "Face or line elements not supported")
            values = []
            if self.num_element_values > 0:
                expect_line(f, "Values:")
                while len(values) < self.num_element_values:
                    line = f.readline()
                    values.extend(map(float, line.split()))
    
            expect_line(f, "Nodes:")
            nodes = []
            while len(nodes) < self.num_nodes:
                line = f.readline()
                nodes.extend(map(int, line.split()))
    
            expect_line(f, "Scale factors:")
            scale_factors = []
            while len(scale_factors) < self.num_scale_factors:
                line = f.readline()
                scale_factors.extend(map(float, line.split()))
    
            self.elements.append(
                    ExelemElement(indices, nodes, values, scale_factors))

    def element_values(self, field_name, component_name, element_num):
        """Return the all field component values at the given element number
        """
        element = self.elements[element_num - 1]
        value_index = 0
        for field in self.fields:
            for component in field.components:
                if component.component_type == 'grid based':
                    component_num_values = np.product(
                            [i + 1 for i in component.divisions])
                    if (field.name == field_name and
                            component.name == str(component_name)):
                        return element.values[value_index:value_index + component_num_values]
                    else:
                        value_index += component_num_values
        raise ValueError("Couldn't find field and component values")


class ExelemField(object):
    """A field definition from an exelem file
    """
    def __init__(self, f, exelem):
        self.exelem = exelem

        # Eg:
        # 1) Geometry, coordinate, rectangular cartesian, #Components=3
        self.components = []

        declaration = f.readline().split(',')
        index, self.name = read_string_regex(f, declaration[0],
                r'([0-9]+)\)\s+([A-Za-z0-9_\s\/]+)')
        self.index = int(index)
        self.num_components = int(read_string_regex(f, declaration[3],
                r'#Components\s*=\s*([0-9]+)$'))
        for i in range(self.num_components):
            component = ExelemComponent(f, self)
            self.components.append(component)

    def __repr__(self):
        return '<ExelemField: %d. "%s", #Components=%d>' % (
                self.index, self.name, self.num_components)


class ExelemComponent(object):
    """A field component definition from an exelem file
    """
    def __init__(self, f, field):
        self.field = field
        # Eg:
        # x. l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.
        #   #Nodes= 8
        declaration = f.readline().strip().split(',')
        self.name = read_string_regex(f, declaration[0], r'^([a-zA-Z0-9]+)\.')
        self.component_type = declaration[2].strip().strip('.')
        if self.component_type == 'standard node based':
            self._read_nodal_component(f)
        elif self.component_type == 'grid based':
            self._read_grid_component(f)
        elif self.component_type == 'data point based':
            self._read_data_component(f)
        else:
            raise ExfileError(f, "Unsupported component type: %s" % 
                    self.component_type)

    def _read_nodal_component(self, f):
        self.num_nodes = int(read_regex(f, r'#Nodes\s*=\s*([0-9]+)'))
        self.node_num_values = {}
        self.value_indices = {}
        self.scale_factor_indices = {}
        for node in range(1, self.num_nodes + 1):
            self.node_num_values[node] = read_regex(f,
                    r'[0-9]+\.\s*#Values\s*=\s*([0-9]+)')
            value_indices_line = f.readline().strip()
            scale_factor_indices_line = f.readline().strip()
            self.value_indices[node] = map(int,
                    value_indices_line.split(':')[1].strip().split())
            self.scale_factor_indices[node] = map(int,
                    scale_factor_indices_line.split(':')[1].strip().split())

    def _read_grid_component(self, f):
        grids = f.readline()
        divisions = [d.strip().split('=')[1] for d in grids.split(',')]
        self.divisions = map(int, divisions)

    def _read_data_component(self, f):
        grids = f.readline()

    def __repr__(self):
        return '<ExnodeComponent, "%s": "%s">' % (
                self.component_type, self.name)


class ExelemElement(object):
    """An element from an exelem file
    """
    def __init__(self, indices, nodes, values, scale_factors):
        self.indices = indices
        self.number = indices[0]
        self.nodes = nodes
        self.values = values
        self.scale_factors = scale_factors

    def __repr__(self):
        return "<ExelemElement %d>" % self.indices[0]

    def __str__(self):
        node_list = " ".join(str(n) for n in self.nodes)
        return "%d: %s" % (self.indices[0], node_list)


class FileWithLineNumber(object):
    def __init__ (self, path, *args):
        if path.endswith('.gz'):
            self.file = gzip.open(path, *args)
        else:
            self.file = open(path, *args)
        self.linenum = 0
        self.prev_pos = self.file.tell()
        self.tprev_pos, self.cur_pos = self.prev_pos, self.prev_pos

    def __enter__ (self):
        return self

    def readline(self):
        self.linenum += 1
        line = self.file.readline()
        self.tprev_pos, self.prev_pos, self.cur_pos = self.prev_pos, self.cur_pos, self.file.tell()
        return line
    
    def rollbacktwice(self):
        self.file.seek(self.tprev_pos)
        self.prev_pos = self.file.tell()
        self.tprev_pos, self.cur_pos = self.prev_pos, self.prev_pos

    def rollback(self):
        self.file.seek(self.prev_pos)
        self.prev_pos = self.file.tell()
        self.tprev_pos, self.cur_pos = self.prev_pos, self.prev_pos

    def __exit__ (self, exc_type, exc_value, traceback):
        self.file.close()


class NodeNotFound(KeyError):
    pass


class ExfileError(ValueError):
    def __init__(self, file, message):
        new_message = "Line %d: %s" % (file.linenum, message)
        super(ExfileError, self).__init__(new_message)


def expect_line(f, expected):
    line = f.readline().strip()
    if line != expected:
        raise ExfileError(f, "Expected '%s', got '%s'" % (expected, line))


def read_regex(f, regex):
    line = f.readline().strip()
    match = re.search(regex, line)
    if match is None:
        raise ExfileError(f, "Expected '%s', got '%s'" % (regex, line))
    if len(match.groups()) == 1:
        return match.group(1)
    else:
        return match.groups()


def read_string_regex(f, string, regex):
    match = re.search(regex, string)
    if match is None:
        raise ExfileError(f, "Expected '%s', got '%s'" % (regex, string))
    if len(match.groups()) == 1:
        return match.group(1)
    else:
        return match.groups()
