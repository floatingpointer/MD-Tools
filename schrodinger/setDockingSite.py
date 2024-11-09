"""
Module/Script to set a "Docking Site" column in Maestro

The maxDistance variable should be set to the threshold outside of which the docking site
will be marked as '0'.

ligandMatches is a dictionary {str:[str,int]}: the key str is reserved for
future use in more complex situations; the value list should contain a SMARTS
identifier str such that each ligand should match to only one SMARTS string
(but the same SMARTS string can be used for multiple ligands).
The second item in the value list is the SMARTS index for the atom used to
calculate distances to the target docking site.

The dockingSiteAtoms() function defines the dictionary objects used that
establish the atoms in the receptor structure that determine the docking sites
(int:str where the int is the atom number and the string is an identifier
that is written to the "Docking Site" column).

The module may be run from within Maestro by loading the module in the GUI's
Python Shell and calling the fromActiveProject() function after selecting the
receptor for the workspace (so it is displayed) and highlight-selecting all
docked ligands (project-table selection) before calling the function.
The appropriate docking site identifier will be written to the "Docking Site"
user-defined column.

Disclaimer: This script has not been tested extensively. Backup your data
before use, and use only at your own risk.
"""

import schrodinger as SD, \
    schrodinger.adapter, \
    schrodinger.maestro, \
    schrodinger.structure, \
    schrodinger.structutils.measure

maxDistance = 6.0 # max distance in angstroms from active site, otherwise '0'

ligandMatches = {
    "U":  ["[H]c1c(=O)n([H])c(=O)nc1[H]", 3] # uracil, C4
}

def dockingSiteAtoms(receptor) -> {int:str}:
    """
    Given a receptor, returns a dictionary of index:name pairs
    for atoms that represent docking sites in the receptor.
    """
    retval = {}
    r = receptor.title.lower()

    if "7o4e" in r:
        retval = {
            2100:   "A:Zn1",
            2101:   "A:Zn2",
            2102:   "B:Zn1",
            2103:   "B:Zn2"
        }
    return retval

def fromActiveProject():
    """
    Use active project to define receptor and ligands. The receptor is defined
    by selection in the workspace (i.e., displayed in the hierarchy) and the
    ligands are selected in the project entry list.
    If more than one receptor is selected, nothing happens (no determination
    of docking site is done, and no property values are updated).
    """
    ligands = receptor = None
    w = SD.maestro.get_included_entries()
    p = [(r, r.getStructure()) for r in SD.maestro.project_table_get().selected_rows]

    receptor = w
    ligands = [l[1] for l in p]
    
    if len(receptor) > 1:
        return
    
    receptor = receptor[0]

    ligands = updateDockingSites(receptor, ligands, dockingSiteAtoms(receptor))
    
    for i in range(len(ligands)):
        p[i][0].setStructure(ligands[i])

    return

def updateDockingSites(receptor, ligands, siteAtoms) -> [SD.structure.Structure]:
    """
    Identifies the docking site by comparing distances to each site.
    Returns the modified list of ligands.
    """

    for l in ligands:
        # SMARTS indexes are 1-based
        for m in ligandMatches:
            smarts = ligandMatches[m][0]
            ligatomSmartsIndex = ligandMatches[m][1]-1

            try:
                al = l.atom[
                    SD.adapter.evaluate_smarts(l,smarts,True)[0][ligatomSmartsIndex]
                ]
            except SD.adapter.SMARTSParseError:
                continue # next ligandMatch

            # matched ligand to SMARTS, so find site by min distance

            site = None
            dmin = 0
            
            for ratom in siteAtoms:
                # ratom is the atom number in the receptor
                ar = receptor.atom[ratom]
                d = SD.structutils.measure.measure_distance(ar, al)
                if site == None or d < dmin:
                    dmin = d
                    site = ratom
            s = '0'
            if dmin < maxDistance:
                s = siteAtoms[site]
            l.property['s_user_Docking_Site'] = s

    return ligands

if __name__ == '__main__':

    import sys
    from pathlib import Path

    infile = Path(sys.argv[1]).absolute()
    outfile = infile
    
    sr = SD.structure.StructureReader(infile)
    
    s = [s for s in sr]

    sr.close()

    receptor = s[0]
    ligands = s[1:]
     
    ligands = updateDockingSites(receptor, ligands, dockingSiteAtoms(receptor))

    sw = SD.structure.StructureWriter(outfile)
    sw.append(receptor)
    sw.extend(ligands)

