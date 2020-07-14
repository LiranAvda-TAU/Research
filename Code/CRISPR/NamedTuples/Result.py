from recordtype import recordtype

Result = recordtype('Result', [('success', False),
                               ('enzymes', None),
                               ('sense_crrnas', []),
                               ('anti_sense_crrnas', []),
                               ('crRNA', None),
                               ('crRNA_strand', 0),
                               ('pam_sites', None),
                               ('ssODN_strand', 0),
                               ('no_extra_inserted_mutations_sites', []),
                               ('inserted_sites', []),
                               ('no_extra_removed_mutations_sites', []),
                               ('removed_sites', [])])
