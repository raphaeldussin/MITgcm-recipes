

def define_faces_connection_ASTE():
    """ define the faces connection for use with xgcm """

    face_connections = {'face':
                    {0: {'Y':  ( None, (1, 'Y', False))},
                     1: {'Y':  ((0, 'Y', False), (2, 'X', False))},
                     2: {'X':  ((1, 'Y', False), (5, 'X', False)),
                         'Y':  (None, (4, 'X', False))},
                     3: {'X':  ((2, 'X', False),  None)},
                     4: {'X':  ((2, 'Y', False),  (5, 'X', False))},
                     5: {'X':  ((4, 'Y', False),  None)} }
    return face_connections
