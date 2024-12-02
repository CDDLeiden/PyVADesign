#!/usr/bin/env python3

import os
import sys
import matplotlib.pyplot as plt
import biotite.sequence.graphics as graphics
from biotite.sequence import Feature, Location, Annotation
from dna_features_viewer import GraphicFeature, GraphicRecord

from .mutation import Mutation
from .sequence import Vector, Gene
from .DNABlocks import DNABlockDesign



class Plot:
    """
    Plotting class for the DNABlocks design.

    Parameters
    ----------
    DNABlocks_design_instance : DNABlockDesign
        DNABlocks design instance.
    mutation_instance : Mutation
        Mutation instance.
    vector_instance : Vector
        Vector instance.
    gene_instance : Gene
        Gene instance.
    output_dir : str
        Output directory.
    show : bool
        Show plot.
    save : bool
        Save plot.
    """
    def __init__(self, 
                 DNABlocks_design_instance: DNABlockDesign,
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 output_dir: str = None,

                 show: bool = False,
                 save: bool = True):
        
        self.DNABlocks_design_instance = DNABlocks_design_instance
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.gene_instance = gene_instance
        self.output_dir = output_dir
        self.show = show
        self.save = save

    def save_plot(self, fig, filename, dpi=100, bbox_inches='tight', transparent=True):
        """Save plot to file."""
        fig.savefig(os.path.join(self.output_dir, filename), dpi=dpi, bbox_inches=bbox_inches, transparent=transparent)

    def plot_histogram_mutations(self, figure_width=8, figure_length=5, filename='histogram_mutations.png'):
        """Show for every DNABlock the number of mutations you can make."""
        counts = self.DNABlocks_design_instance.count_mutations_per_DNABlock()
        fig, ax = plt.subplots(figsize=(figure_width, figure_length))
        labels = list(counts.keys())
        colors = list(self.DNABlocks_design_instance.DNABlock_colors.values())[:len(counts)]
        ax.bar(range(len(counts)), list(counts.values()), align='center', color=colors)
        plt.xticks(range(len(counts)), labels, rotation=90)
        ax.set_ylabel('Number of mutants per DNABlock')
        ax.set_title(f'Number of mutants per DNABlock')
        if self.save:
            self.save_plot(fig, filename)
        if self.show:
            plt.show()
        else:
            plt.close()
        
    def plot_mutation_legend(self, legend_alpha=0.2, font_size='x-large', marker_size=10, linestyle='None', marker='o', loc='center', bbox_to_anchor=(0.5, 0.5)):
        """Legend for DNABlocks plot."""
        # Create an empty plot with no data points
        fig, ax = plt.subplots(figsize=(3, 2))
        # Add mutation type colors to the legend
        handles = []
        for k, v in self.mutation_instance.colors.items():
            handle = plt.Line2D([0], [0], marker=marker, color=f'{v}', label=f'{k}', markersize=marker_size, linestyle=linestyle, alpha=legend_alpha)
            handles.append(handle) 
        legend = ax.legend(handles=handles, loc=loc, bbox_to_anchor=bbox_to_anchor, fontsize=font_size, framealpha=legend_alpha)
        # Hide the axes
        ax.axis('off')
        if self.show:
            plt.show()
        else:
            plt.close()
        if self.save:
            self.save_plot(fig, 'legend.png')

    def plot_DNABlocks_mutations(self, 
                               plot_DNABlocks=True, 
                               plot_mutations=True, 
                               figure_width=20, 
                               figure_length=10):
        """Plot mutations and selected DNABlocks."""
        # Check if mutations and sequence are available
        if self.mutation_instance.mutations is None:
            raise ValueError("No mutations found. Please run the mutation class first.")
        if self.gene_instance.sequence is None:
            raise ValueError("No sequence found. Please run the sequence class first.")
        features = []
        # Add GoI to plot
        features.append(GraphicFeature(start=0, 
                                       end=len(self.gene_instance.sequence) +3, 
                                       strand=+1, 
                                       color=Vector(Gene).color, 
                                       label=f"{self.gene_instance.seqid}"))
        # Add mutations to plot
        if plot_mutations:
            for i in self.mutation_instance.silent_mutations:
                features.append(GraphicFeature(start=int(i.idx_dna[0]), 
                                               end=int(i.idx_dna[0]) + 3,
                                               strand=+1, 
                                               color=self.mutation_instance.colors[i.type], 
                                               label=f"{i.name}"))

            for num, mut in enumerate(self.mutation_instance.mutations):
                if not mut.is_multiplemutation:
                    features.append(GraphicFeature(start=int(mut.idx_dna[0]), 
                                                   end=int(mut.idx_dna[0]) + 3,
                                                   strand=+1, 
                                                   color=self.mutation_instance.colors[mut.type],
                                                   label=f"{''.join(mut.mutation)}"))

                else:
                    for num, _ in enumerate(mut.idx_dna):
                        features.append(GraphicFeature(start=int(mut.idx_dna[num]), 
                                                       end=int(mut.idx_dna[num]) + 3,
                                                       strand=+1, 
                                                       color=self.mutation_instance.colors['Combined'], 
                                                       label=f"{mut.mutation[num]}"))
                        
        # Add DNABlocks to plot
        if plot_DNABlocks:
            for num, i in enumerate(self.DNABlocks_design_instance.wt_DNABlocks):

                features.append(GraphicFeature(start=i.bin_start, 
                                               end=i.bin_end,
                                               strand=+1, 
                                               color=self.DNABlocks_design_instance.DNABlock_colors[num + 1],
                                               label=f"{i.name}"))

            
        record = GraphicRecord(sequence_length=len(self.gene_instance.sequence), features=features)
        fig_size = (figure_length, figure_width)
        fig, ax = plt.subplots(figsize=fig_size) 
        record.plot(ax=ax, figure_width=figure_width)
        if self.show:
            plt.show()
        else:
            plt.close()

        if plot_DNABlocks and plot_mutations and self.save:
            if self.DNABlocks_design_instance.cost_optimization:
                self.save_plot(fig, f'DNABlocks_{self.gene_instance.seqid}_N{self.mutation_instance.n_mutants}_cost.png')
            elif self.DNABlocks_design_instance.amount_optimization:
                self.save_plot(fig, f'DNABlocks_{self.gene_instance.seqid}_N{self.mutation_instance.n_mutants}_amount.png')
        if not plot_DNABlocks:
            self.save_plot(fig, f'mutations_{self.gene_instance.seqid}_N{self.mutation_instance.n_mutants}_{self.DNABlocks_design_instance.optimization_method}.png')

    def extract_snapgene_features(self):
        """Extract features from a snapgene vector file and return them as a list of biotite sequence features."""
        if self.vector_instance.vector is None:
            raise ValueError("No vector found. Please run the sequence class first.")

        vector_features = []
        for i in self.vector_instance.vector.features:
            feature_type = i.type
            quals = i.qualifiers
            label_qualifier = i.qualifiers.get('label', '')
            locations = [Location(int(i.location.start), int(i.location.end))]

            if label_qualifier:
                qual = {"label": ' '.join(label_qualifier)}
            else:
                qual = {}

            if (i.type == "rep_origin") or (i.type == "protein_bind") or (i.type == "terminator"):
                if quals['note']:
                    qual = {"note": ' '.join(quals['label'])}

            elif (i.type == "CDS") or (i.type == "gene"):
                if quals['label']:
                    qual = {"product": ' '.join(quals['label'])}
            elif (i.type == "regulatory"):
                if quals['label']:
                    qual = {"label": ' '.join(quals['label'])}
            elif (i.type == "misc_feature"):
                if quals['label']:
                    qual = {"note": ' '.join(quals['label'])}

            vector_features.append(Feature(feature_type, locations, qual))
        # Add organism source
        if self.vector_instance.organism:
                vector_features.append(Feature("source", [Location(0, len(self.vector_instance.vector.seq))], {"organism": f"{self.vector_instance.organism}"}))
        return vector_features

    def plot_vector(self, figsize=(8,8), fontsize=10):
        """Plot a plasmid map of the vector sequence."""
        features = self.extract_snapgene_features()
        annotation = Annotation(features)
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="polar")
        graphics.plot_plasmid_map(
            ax, 
            annotation, 
            plasmid_size=len(self.vector_instance.vector.seq), 
            label=f"{self.vector_instance.vector.name}",
            label_properties={"fontsize": fontsize})
        ticks = ax.get_xticks()
        # Only show a few ticks
        ax.set_xticks(ticks[::2])
        fig.tight_layout()
        if self.save:
            self.save_plot(fig, f"plasmid_map.png")
        if self.show:
            plt.show()
        else:
            plt.close()
        return ax, fig