newick = ";"

nexml = """<nexml version="0.9" xmlns="http://www.nexml.org/2009">
  <trees id="emptyTreeSet">
    <tree id="emptyTree"/>
  </trees>
</nexml>"""

nexus = """#NEXUS
BEGIN TREES;
  TREE emptyTree = ();
END;"""

lookup = {
    "newick": newick,
    "nexml": nexml,
    "nexus": nexus,
}
