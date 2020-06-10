from recordtype import recordtype

Result = recordtype('Result', [('success', False),
                               ('crRNA', None),
                               ('ssODN_strand', 0),
                               ('no_extra_inserted_mutations_sites', []),
                               ('inserted_sites', []),
                               ('inserted_favourite', False),
                               ('no_extra_removed_mutations_sites', []),
                               ('removed_sites', []),
                               ('removed_favourite', False)])
