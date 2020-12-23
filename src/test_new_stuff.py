#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 12:12:13 2020

@author: elia
"""


import json

def read_json_file(filename: str) -> dict:
    """Reads a JSON file and returns a dictionary with the content

    Args:
        filename: Name of the json file to read

    Returns:
        Dictionary with the json file info

    """

    with open(filename) as json_content:

        return json.load(json_content)


JSON_CONFIG_FILENAME = "config.json"
config = read_json_file(JSON_CONFIG_FILENAME)

configOrganism = config["organism"]
configOrganismFactory = config["organismFactory"]
configConnector = config["connector"]
configPssm = config["pssm"]



# ---------------------------------TEST----------------------------------------

organism_factory = OrganismFactory(
        configOrganism, configOrganismFactory, configConnector, configPssm
    )


# Test get_organism
parent_1 = organism_factory.get_organism()  # Parent 1

# Test 'count' methods
print('Number of recognizers:', parent_1.count_recognizers())
print('Number of recognizers:', len(parent_1.recognizers))

print('Number of connectors:', parent_1.count_connectors())
print('Number of connectors:', len(parent_1.connectors))


# Create another organism
parent_2 = organism_factory.get_organism()  # Parent 2
print('Number of recognizers:', parent_2.count_recognizers())
print('Number of connectors:', parent_2.count_connectors())


# Test 'print' method
print("\n")
parent_1.print()  # Print parent_1
print("\n")
parent_2.print()  # Print parent_2


# Recombination process
# Test get_children
child_1, child_2 = organism_factory.get_children(parent_1, parent_2)

# Print children
print("\n")
child_1.print()  # Print child_1
print("\n")
child_2.print()  # Print child_2





# Test 'export' method
parent_1.export("test_export_org.txt")
child_1.export("test_export_org.txt")


















