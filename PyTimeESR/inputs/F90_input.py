
import re

from .default_inputs import *

class _complex_format(str):
    """Custom string class to format complex numbers.
    """

    def format(self):
        return '({:.8f},{:.8f})'.format(self.real, self.imag)              

class F90Input(): 
    """Parent class to write fortran input.
    """
    
    fmt = {
        tfloat: '{:.8f} ',
        tcomplex: _complex_format,
        tint: '{:d} ',
        tstr: '{} ',
        tbool: '{} ',
    }
    
    def __init__(self, line_lenght, padding_lenght):
        self.line_lenght = line_lenght
        self.padding_lenght = padding_lenght

    @staticmethod
    def check_dictionary(dictionary, required_keys, dict_label):
        """Check if the dictionary has the required keys and types.
        Args
        -----
        dictionary: dict
            Dictionary to check.
        required_keys: dict
            Dictionary with the required keys and types. if the type is a list, the first element
            is the type of the list, the second element type of elememts,
            and the third element is the number of elements in the list.
        dict_label: str
            Label of the dictionary to check. to be used in the error message.
        """
        assert isinstance(dictionary, dict), "Hamiltonian input should be a dictionary"
        for key, value in required_keys.items(): 
            assert key in dictionary, f"Missing key {key} in {dict_label}"

            if not isinstance(value, tlist):
                # if not a list check if the type is correct
                assert isinstance(dictionary[key], value), f"Key {key} should be of type {value} in {dict_label}"
                continue
            
            # if itterable check if the type is correct
            assert isinstance(dictionary[key], value[0]), f"{key} should be of type {value} in {dict_label}"
            # check if the elements of the list are of the correct type
            assert all(isinstance(i, value[1]) for i in dictionary[key]), f"Elements of {key} should be of type {value} in {dict_label}"
            if value[2] is None:
                continue
            # check if the number of elements is correct if it is not fixed
            assert len(dictionary[key]) == value[2], f"Key {key} should have {value[2]} elements in {dict_label}"

    def create_header(self, header: str, seperator: str):
        """Create the header for the input file.
        Args
        -----
        header: str
            Header to add to the input file.
        seperator: str
            Seperator to add to the input file.
        lenghth: int
            Length of the header. Default is 80.
        """

        l = int((self.line_lenght - len(header))/2)
        line = seperator * l + header + seperator * l + '\n'
        return line
    
    @staticmethod
    def bool2string(value: bool):
        string = '.true.' if value else '.false.'
        return string
    
    @staticmethod
    def string2bool(value: str):
        if value == '.true.':
            return True
        elif value == '.false.':
            return False
        else:
            raise ValueError(f"Invalid boolean string: {value}. Use '.true.' or '.false.'.")
        
    @staticmethod
    def load_complex(line: str) -> complex:
        """Load a complex number from a line of the input file.
        format: (Re, Im)

        Args
        -----
        line: str
            Line from the input file containing a complex number.

        Returns
        -------
        complex: 
            Complex number loaded from the line.
        """
        string = re.split(r'[()\s,]', line)[1:3]
        z = complex(*tuple(map(float, string)))
        return z

    def input_line(self, value, comment: str = ''):
        """Create a line for the input file.
        Args
        -----
        value: str
            Value to add to the input file.
        comment: str
            Comment to add to the input file.
        """
        # turn everything into a list
        if not isinstance(value,tlist):
            value = [value]
        
        
        # convert bools to strings
        value = [self.bool2string(val) if isinstance(val, tbool) else val for val in value]
        t = type(value[0])
    
        # convert t the specified type into a generic type ttype
        ttype = tfloat if (t in tfloat) else tint if (t in tint) else tstr if (t in tstr) else tcomplex if (t in tcomplex) else None
        assert ttype is not None, f"Type {t} not supported in input line"
    
        line = ' '.join([self.fmt[ttype].format(i) for i in value])
        pad = self.padding_lenght - len(line)
        pad = 0 if pad < 0 else pad
        line += ' ' * pad + '! ' + comment + '\n'
        return line
    
    def load_output(self, output_dict: str):
        """Load the output from a file.
        
        Args
        -----
        output_dict: str
            Output file to load.
        
        Returns
        -------
        dict: Dictionary with the output.
        """
        
        return {}
