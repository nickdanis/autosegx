# autosegx: Autosegmental representations using networkx

The `autosegx` module provides some basic tools for working with autosegmental representations (for phonology, and beyond) building off the `networkx` package for python. Autosegmental structure are defined as graphs, with connections between nodes, and then these structures can be compared, searched for natural classes, evaluated against constraints, and more. It started as a basic tool to find the predicted natural classes for a specific theory of representation and to compare these across theories, and is slowly turning into a more general set of tools.

## Basic usage