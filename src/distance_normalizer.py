import numpy as np

print('[INFO] Importing pyrosetta')
import pyrosetta

print('[INFO] Initializing rosetta')
pyrosetta.init('-out:level 0 -ignore_unrecognized_res')


def apply_distances(distances, abmodel):
    print('[INFO] Normalizing with NATIVE distances')
    diff = np.array([0.0, 0.0, 0.0])

    new_model = []

    for index, coord in enumerate(abmodel[:-1]):
        native_distance = distances[index]
        p1 = np.array(coord[1:])
        p2 = np.array(abmodel[index + 1][1:])

        pos_diff = p2 - p1
        offset = pos_diff * (native_distance - 1)
        diff += offset
        p2 += diff

        new_model.append([coord[0], p2[0], p2[1], p2[2]])

    p2 = np.array(abmodel[-1][1:]) + diff

    new_model.append([abmodel[-1][0], p2[0], p2[1], p2[2]])

    return new_model


def dump_abmodel_as_pdb(target, abmodel):
    name = '%s_3dab.pdb' % target
    print('[INFO] Dumping adjusted AB model to %s' % name)

    with open(name, 'wt') as f:
        for index, data in enumerate(abmodel):
            amino_code, pos_x, pos_y, pos_z = data

            f.write("ATOM    %3d  CA  %s A %3d      %6.3f  %6.3f  %6.3f  1.00  0.00           C  \n" % (
                index + 1, amino_code, index + 1, pos_x, pos_y, pos_z
            ))

        f.write('TER\n')


def get_distances(pose, abmodel):
    print('[INFO] Measuring distances')

    positions = []
    distances = []

    for i in range(len(abmodel)):
        coords = pose.residue(i + 1).xyz('CA')
        positions.append((coords.x, coords.y, coords.z))

    for k, p1 in enumerate(positions[:-1]):
        p1 = np.array(p1)
        p2 = np.array(positions[k + 1])

        dist = np.linalg.norm(p2 - p1)
        distances.append(dist)

    return distances


def normalize_distances(native_filename, abmodel):
    native_pose = pyrosetta.pose_from_pdb(native_filename)
    distances = get_distances(native_pose, abmodel)
    normalized = apply_distances(distances, abmodel)
    dump_abmodel_as_pdb(native_filename.split('.')[0], normalized)

