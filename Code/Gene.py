class Gene:

    def __init__(self, name, species):
        self.name = name
        self.species = species
        self.id = None
        self.ncbi_id = None

    def set_id(self, gene_id):
        self.id = gene_id

    def set_ncbi_id(self, ncbi_id):
        self.ncbi_id = ncbi_id
