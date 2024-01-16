class Eblocks:
    def __init__(self):
        self.eblock_parameters = {"common_param": "shared_value"}
        self.block_sequences = []

    def set_parameters(self, parameters):
        self.eblock_parameters = parameters

    def set_block_sequences(self, block_sequences):
        self.block_sequences = block_sequences

    def display_parameters(self):
        print("Eblocks Parameters:", self.eblock_parameters)
        print("Block Sequences:", self.block_sequences)


class Design:
    def __init__(self, eblocks_instance):
        self.eblocks_instance = eblocks_instance
        self.block_sequences = []

    def design_eblocks(self):
        # Your design logic to generate block_sequences
        self.block_sequences = ["sequence1", "sequence2", "sequence3"]

        # Set the block_sequences in the Eblocks instance
        self.eblocks_instance.set_block_sequences(self.block_sequences)


# Example usage
eblocks_instance = Eblocks()
design_instance = Design(eblocks_instance)

# Design eblocks and set block_sequences in Eblocks instance
design_instance.design_eblocks()

# Display parameters from both Eblocks and Design instances
eblocks_instance.display_parameters()