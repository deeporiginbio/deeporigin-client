# TODO: use openbabel function running on knative to support this
# def convert_file(
#     source_type, source, destination_type, destination=None, add_hydrogens=False
# ):
#     from openbabel import pybel

#     mol_pb_gen = pybel.readfile(str(source_type), source)

#     if not destination:
#         destination = tempfile.mktemp()

#     out_file = pybel.Outputfile(str(destination_type), destination, overwrite=True)

#     for mol_pb in mol_pb_gen:
#         if add_hydrogens:
#             mol_pb.addh()
#         out_file.write(mol_pb)

#     out_file.close()
#     return destination


# def convert_block(source_type, source, destination_type, add_hydrogens=False):
#     source_file = tempfile.mktemp()
#     with open(source_file, "w") as f:
#         f.write(source)

#     destination_file = convert_file(
#         source_type, source_file, destination_type, add_hydrogens=add_hydrogens
#     )

#     with open(destination_file, "r") as f:
#         return f.read()
