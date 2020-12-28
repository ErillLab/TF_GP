#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 12:12:13 2020

@author: elia
"""


from objects.organism_factory import OrganismFactory

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

#------------------------------------------------------------------------------
# Test get_organism
#------------------------------------------------------------------------------

parent_1 = organism_factory.get_organism()  # Parent 1

# Create another organism
parent_2 = organism_factory.get_organism()  # Parent 2


#------------------------------------------------------------------------------
# Test 'print' method
#------------------------------------------------------------------------------

print("\n")
parent_1.print()  # Print parent_1
print("\n")
parent_2.print()  # Print parent_2


#------------------------------------------------------------------------------
# Test 'count' methods
#------------------------------------------------------------------------------

# (parent_1)
# Number of recognizers
print('Number of recognizers:', parent_1.count_recognizers())
# Same thing as
print('Number of recognizers:', len(parent_1.recognizers))

# Number of connectors
print('Number of connectors:', parent_1.count_connectors())
# Same thing as
print('Number of connectors:', len(parent_1.connectors))


# (parent_2)
print('Number of recognizers:', parent_2.count_recognizers())
print('Number of connectors:', parent_2.count_connectors())


#------------------------------------------------------------------------------
# Test get_children
#------------------------------------------------------------------------------

# Recombination process
child_1, child_2 = organism_factory.get_children(parent_1, parent_2)

# Print children
print("\n")
child_1.print()  # Print child_1
print("\n")
child_2.print()  # Print child_2


#------------------------------------------------------------------------------
# Test 'export' method
#------------------------------------------------------------------------------

parent_1.export("test_export_parent1.txt")
child_1.export("test_export_child1.txt")

#------------------------------------------------------------------------------
# Test json import/export
#------------------------------------------------------------------------------

# Three random organisms
org_1 = organism_factory.get_organism()
org_2 = organism_factory.get_organism()
org_3 = organism_factory.get_organism()

print("\n")
org_1.print()  # Print org_1
print("\n")
org_2.print()  # Print org_2
print("\n")
org_3.print()  # Print org_3

# Export organisms
organism_factory.export_organisms([org_1, org_2, org_3], "three_random_org_json_test.json")
# Import organisms
org_list = organism_factory.import_organisms("three_random_org_json_test.json")

print("\n")
org_list[0].print()  # Print org_1
print("\n")
org_list[1].print()  # Print org_2
print("\n")
org_list[2].print()  # Print org_3


#------------------------------------------------------------------------------
# Test PSSM scores
#------------------------------------------------------------------------------

new_pssm = organism_factory.create_pssm(4)  # PSSM

new_pssm.print()

new_pssm.get_score("ccag")
new_pssm.get_score("ctgg")








