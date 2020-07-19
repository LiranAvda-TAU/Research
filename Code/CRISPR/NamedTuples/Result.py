from recordtype import recordtype

Result = recordtype('Result', [('query_data', None),
                               ('success', False),
                               ('enzymes', None),
                               ('sense_crrnas', []),
                               ('anti_sense_crrnas', []),
                               ('crRNA', None),
                               ('crRNA_strand', 0),
                               ('pam_sites', None),
                               ('ssODN_strand', 0),
                               ('no_extra_inserted_mutations', []),
                               ('inserted_mutations', []),
                               ('no_extra_removed_mutations', []),
                               ('removed_mutations', [])])
