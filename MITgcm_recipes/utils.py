

def define_faces_connection_ASTE():
    """ define the faces connection for use with xgcm """
    face_connections = {'face':
                    {0: {'X':  ((5, 'Y', False), None),
                         'Y':  (None,             (1, 'Y', False))},
                     1: {'X':  ((4, 'Y', False), None),
                         'Y':  ((0, 'Y', False),  (2, 'X', False))},
                     2: {'X':  ((1, 'Y', False),  (3, 'X', False)),
                         'Y':  (None,             (4, 'X', False))},
                     3: {'X':  ((2, 'X', False), None),
                         'Y':  (None,            None)},
                     4: {'X':  ((2, 'Y', False),  (5, 'X', False)),
                         'Y':  (None,             (1, 'X', False))},
                     5: {'X':  ((4, 'X', False), None),
                         'Y':  (None,             (0, 'X', False))}}}
    return face_connections
