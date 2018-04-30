# -*- coding: utf-8 -*-
"""
Classes :
* Coord
* ProfilTravers
* SuiteProfilsTravers
* LignesContrainte
* Lit
* MeshConstructor
"""

# Compatibility with Python2
from __future__ import print_function

from copy import deepcopy
from math import ceil
import numpy as np
from numpy.lib.recfunctions import append_fields
import pandas as pd
from shapely.geometry import LineString, Point
import sys

from geom.dataset import BlueKenueRead_i2s, BlueKenueRead_i3s, BlueKenueWrite_i2s


# FUNCTIONS
def get_intersections(linestrings):
    """
    Chercher les intersections entre toutes les lignes linestrings
    linestrings <list <LineString>>
    """
    intersections = []
    for i, l1 in enumerate(linestrings):
        for j, l2 in enumerate(linestrings[i+1:]):
            if l1.intersects(l2):
                intersections.append((i, j+i+1))
    return intersections


def float_vars(varname_list):
    return [(label, np.float) for label in varname_list]


def strictly_increasing(L):
    """Check if array L is sorted"""
    return all(x < y for x, y in zip(L, L[1:]))


def get_axe_hydraulique(infile_axe):
    """Extraction de l'unique ligne du fichier i2s d'entrée"""
    with BlueKenueRead_i2s(infile_axe) as in_i2s:
        in_i2s.read_header()
        lines = [line for value, line in in_i2s.iter_on_polylines()]
    nb_lines = len(lines)
    if nb_lines != 1:
        sys.exit("Le fichier '{}' contient {} polylignes au lieu d'une seule pour définir l'axe hydraulique".format(infile_axe, nb_lines))
    return lines[0]


# CLASSES
class Droite:
    """
    Droite: description paramétrique d'une droite
    """
    def __init__(self, x_P1, y_P1, x_P2, y_P2):
        """Équation paramétrique : a*x + b*y = 1"""
        self.P1 = (x_P1, y_P1)
        self.P2 = (x_P2, y_P2)
        self.a = (y_P2 - y_P1)/(y_P2*x_P1 - y_P1*x_P2)
        self.b = (x_P2 - x_P1)/(x_P2*y_P1 - x_P1*y_P2)

    def intersection(self, other):
        """
        Recherche le point d'intersection entre self et other
        /!\ Le script devrait planter si les droites sont parrallèles...
        """
        yC = (self.a - other.a)/(self.a*other.b - other.a*self.b)
        xC = (1 - self.b*yC)/self.a
        return xC, yC

    def ratio_along_line(self, xM, yM):
        """
        Calcul de la position relative du point M par rapport à la droite initiale
        (points P1 et P2 utilisée pour construire la droite)
        0 (M=P1) <= ratio <= 1 (M=P2)
        """
        if self.P2[0] != self.P1[0]:
            ratio = (xM - self.P1[0])/(self.P2[0] - self.P1[0])
        else:
            ratio = (yM - self.P1[1])/(self.P2[1] - self.P1[1])
        return max(min(ratio, 1.), 0.)


class Coord:
    """
    Coord: ensemble de points 2D/3D consécutifs (ex: polylignes)

    ### Attributs
    - array

    ### Méthodes
    - compute_Xt
    - compute_xt
    - compute_xp
    - move_between_targets
    """
    LABELS = ['X', 'Y', 'Z', 'Xt', 'xt']  # nevers used...

    def __init__(self, array, vars2add, remove_duplicates=False):
        """
        array: X, Y (Z is optional)
        Add Xt and xt columns if not already present
        """
        self.array = array

        vars = list(self.array.dtype.fields)
        if 'X' not in vars and 'Y' not in vars:
            sys.exit("Columns X and Y are compulsary. Following columns were found: {}".format(vars))

        if 'Xt' in vars2add:
            self.compute_Xt()
        if 'xt' in vars2add:
            self.compute_xt()

        if remove_duplicates:
            if not strictly_increasing(self.array['Xt']):
                print("ATTENTION : Des points doublons sont éliminés")
                # Suppression des doublons (points superposés dans la polyligne)
                points_a_conserver = np.ediff1d(self.array['Xt'], to_begin=1.) != 0.
                self.array = self.array[points_a_conserver]

    def compute_Xt(self):
        """
        Calcul de l'abscisse curviligne (distance 2D cumulée)
        """
        Xt = np.sqrt(np.power(np.ediff1d(self.array['X'], to_begin=0.), 2) +
                     np.power(np.ediff1d(self.array['Y'], to_begin=0.), 2))
        Xt = Xt.cumsum()
        self.array = append_fields(self.array, 'Xt', Xt, usemask=False)

    def compute_xt(self):
        """
        Calcul de l'asbcisse curviligne adimensionnée (de 0 à 1)
        /!\ La colonne 'Xt' doit déjà exister
        """
        if len(self.array) > 1:
            xt = (self.array['Xt'] - self.array['Xt'][0])/(self.array['Xt'][-1] - self.array['Xt'][0])
        else:
            xt = np.empty(len(self.array))
            xt.fill(-999.)
        self.array = append_fields(self.array, 'xt', xt, usemask=False)

    def move_between_targets(self, p1, p2):
        """
        Déplace la ligne pour la faire correspondre aux limites voulues
        (par translation pondérée linéairement de p1 à p2)
        p1, p2 <Point>: point 1 (origine ligne) et 2 (fin ligne)
        /!\ L'objet est modifié directement et la méthode ne retourne rien
        """
        array = deepcopy(self.array)
        xt = array['xt']
        self.array['X'] = array['X'] + (1-xt)*(p1.x-array['X'][0]) + xt*(p2.x-array['X'][-1])
        self.array['Y'] = array['Y'] + (1-xt)*(p1.y-array['Y'][0]) + xt*(p2.y-array['Y'][-1])

    def compute_xp(self):
        """Calcul de la distance projetée adimensionnée sur droite début->fin"""
        self.trace = LineString([self.array[['X','Y']][0], self.array[['X','Y']][-1]])
        Xp = np.empty(len(self.array), dtype=np.float)

        for i, row in enumerate(self.array):
            point = Point(list(row)[:2])
            Xp[i] = self.trace.project(point)

        if not strictly_increasing(Xp):
            sys.exit("L'abscisse projetée n'est pas strictement croissante le long du profil. Le profil n'est pas suffisamment droit...")

        xp = Xp/Xp[-1]
        self.array = append_fields(self.array, 'xp', xp, usemask=False)

    def convert_as_array(self):
        return np.column_stack((self.array['X'], self.array['Y'], self.array['Z']))

    def convert_as_linestring(self):
        return LineString(self.convert_as_array())


class ProfilTravers:
    """
    ProfilTravers: représente un profil en travers

    ### Attributs
    - id <integer>: identifiant unique (numérotation automatique commençant à 0)
    - coord <Coord>
    - geom: objet géometrique LineString
    - limites: dictionnaire du type: {id_ligne: (Xt_profil, Xt_ligne, intersection.z)}
    - dist_proj_axe (créée par une méthode de <SuiteProfilsTravers>

    ### Méthodes
    - _add_limit
    - get_limit_by_id
    - find_and_add_limit
    - common_limits
    - extraire_lit
    - sort_limites
    - export_profil_csv
    - export_limites_csv
    - calcul_nb_pts_inter
    """
    LABELS = ['X', 'Y', 'Z', 'Xt']

    def __init__(self, id, coord, label):
        """
        Créer un profil à partir d'un semis de points ordonnés (coordonnées X,Y,Z)
        Aucune limite de lit n'est déclarée lorsque l'objet est créé
        """
        self.id = id
        self.label = label
        self.coord = Coord(np.array(coord, dtype=float_vars(['X', 'Y', 'Z'])), ['Xt'])
        self.nb_points = len(self.coord.array)
        self.geom = LineString(coord)  # FIXME: contient toujours les points doublons

        self.limites = pd.DataFrame(columns=['Xt_profil', 'Xt_ligne', 'X', 'Y', 'Z', 'id_pt'])
        self.limites['id_pt'] = self.limites['id_pt'].astype(int)  # FIXME change dtype

    def __repr__(self):
        return "{} #{} ({} points)".format(self.label, self.id, self.nb_points)

    def _add_limit(self, id_ligne, Xt_profil, Xt_ligne, point):
        """Ajoute une nouvelle limite au profil"""
        corr_bug = self.geom.interpolate(Xt_profil).z  # DEBUG: l'altitude du point semble être buggée (ie point.z) pour de grosses données? Il vaut mieux la recalculer avec Xt_profil
        row = pd.DataFrame({'Xt_profil': Xt_profil, 'Xt_ligne': Xt_ligne,
                            'X': point.x, 'Y': point.y, 'Z': corr_bug},
                           index=[id_ligne])
        self.limites = self.limites.append(row)

    def get_limit_by_id(self, id_ligne):
        return self.limites.loc[id_ligne]

    def get_Xt_lignes(self, id1, id2):
        return tuple(self.limites.loc[[id1,id2],'Xt_ligne'])

    def find_and_add_limit(self, ligne_contrainte, dist_max):
        if self.geom.intersects(ligne_contrainte.geom):
            intersection = self.geom.intersection(ligne_contrainte.geom)

            if isinstance(intersection, Point):
                # Calcul des projections
                Xt_profil = self.geom.project(intersection)
                Xt_ligne = ligne_contrainte.geom.project(intersection)
                self._add_limit(ligne_contrainte.id, Xt_profil, Xt_ligne, intersection)
            else:
                sys.exit("L'intersection entre '{}' et '{}' ne correspond pas à un point unique (type = {})".format(self, ligne_contrainte, type(intersection)))

        else:
            distance = self.geom.distance(ligne_contrainte.geom)
            if dist_max is not None:
                if distance < dist_max:
                    # Test de proximité avec tous les points de la ligne de contrainte
                    i = 0
                    for coord in ligne_contrainte.coord:
                        i += 1
                        point = Point(coord)
                        dist = self.geom.distance(point)
                        if dist < dist_max:
                            # On a trouvé le bon point et on ajoute la limite
                            Xt_ligne = ligne_contrainte.geom.project(point)
                            Xt_profil = self.geom.project(point)
                            intersection = self.geom.interpolate(Xt_profil)
                            self._add_limit(ligne_contrainte.id, Xt_profil, Xt_ligne, intersection)
                            print("Ajout de la limite avec la ligne {} après {} itérations (distance = {})".format(ligne_contrainte.id, i, dist))
                            break

    def common_limits(self, other):
        """
        Liste des limites communes entre les profils self et other
        self et other sont de type: <ProfilsTravers> et sont commutables
        """
        return self.limites.index[np.in1d(self.limites.index,
                                          other.limites.index,
                                          assume_unique=True)]

    def extraire_lit(self, id_lit_1, id_lit_2):
        """
        Extraire les coordonnées d'une partie du profil :
        Les points du profil compris entre les deux lits et avec éventuellement les points de bord interpolés
        Retourne un <Lit> avec les colonnes ('X', 'Y', 'Z', 'Xt', 'xt')
        /!\ id_lit_1 et id_lit_2 doivent être dans l'ordre Xt croissant (plante sinon)
        """
        limit1 = self.get_limit_by_id(id_lit_1)
        limit2 = self.get_limit_by_id(id_lit_2)

        Xt1 = limit1['Xt_profil']
        Xt2 = limit2['Xt_profil']

        # Vérifie que les Xt sont croissants du id_lit_1 au id_lit_2
        if Xt1 > Xt2:
            sys.exit("L'ordre des lits ({}, {}) n'est pas par Xt croissant pour le {}".format(id_lit_1, id_lit_2, self))

        Xt_profil = self.coord.array['Xt']
        sub_coord = self.coord.array[np.logical_and(Xt_profil >= Xt1, Xt_profil <= Xt2)]

        # Ajoute le premier point si nécessaire
        if Xt1 not in Xt_profil:
        # if min(abs((Xt_profil - Xt1)) > 0.1):
            row = np.array([(limit1['X'], limit1['Y'], limit1['Z'], Xt1)], dtype=float_vars(ProfilTravers.LABELS))
            sub_coord = np.insert(sub_coord, 0, row)

        # Ajoute le dernier point si nécessaire
        if Xt2 not in Xt_profil:
        # if min(abs((Xt_profil - Xt2)) > 0.1):
            row = np.array([(limit2['X'], limit2['Y'], limit2['Z'], Xt2)], dtype=float_vars(ProfilTravers.LABELS))
            sub_coord = np.append(sub_coord, row)

        # Vérification de l'ordre des points
        if not strictly_increasing(sub_coord['Xt']):
            print("/!\ Les Xt ne sont pas strictement croissants") #FIXME: should not appear
            print(sub_coord['Xt'])
            print("Veuillez vérifier ci-dessus, avec les limites suivantes :")
            print(limit1)
            print(limit2)
            points_a_conserver = np.ediff1d(sub_coord['Xt'], to_begin=1.)!=0.
            sub_coord = sub_coord[points_a_conserver]

        return Lit(sub_coord, ['xt'])

    def extraire_lit_v2(self, id_lit_1, id_lit_2):
        global sub_coord, Xt_profil
        """
        Idem extraire_lit_v2 mais avec le numéro des points
        """
        limit1 = self.get_limit_by_id(id_lit_1)
        limit2 = self.get_limit_by_id(id_lit_2)

        Xt1 = limit1['Xt_profil']
        Xt2 = limit2['Xt_profil']

        # Vérifie que les Xt sont croissants du id_lit_1 au id_lit_2
        if Xt1 > Xt2:
            sys.exit("L'ordre des lits ({}, {}) n'est pas par Xt croissant pour le {}".format(id_lit_1, id_lit_2, self))

        Xt_profil = self.coord.array['Xt']
        rows2keep = np.logical_and(Xt_profil >= Xt1, Xt_profil <= Xt2)
        id_pts = np.nonzero(rows2keep)
        sub_coord = self.coord.array[rows2keep]

        # Ajoute le premier point si nécessaire
        if Xt1 not in Xt_profil:
            row = np.array([(limit1['X'], limit1['Y'], limit1['Z'], Xt1)], dtype=float_vars(ProfilTravers.LABELS))
            sub_coord = np.insert(sub_coord, 0, row)

        # Ajoute le dernier point si nécessaire
        if Xt2 not in Xt_profil:
            row = np.array([(limit2['X'], limit2['Y'], limit2['Z'], Xt2)], dtype=float_vars(ProfilTravers.LABELS))
            sub_coord = np.append(sub_coord, row)

        # Vérification de l'ordre des points
        if not strictly_increasing(sub_coord['Xt']):
            print("/!\ Les Xt ne sont pas strictement croissants") #FIXME: should not appear
            print(sub_coord['Xt'])
            print("Veuillez vérifier ci-dessus, avec les limites suivantes :")
            print(limit1)
            print(limit2)
            raise Exception("FIXME")
            points_a_conserver = np.ediff1d(sub_coord['Xt'], to_begin=1.)!=0.
            sub_coord = sub_coord[points_a_conserver]

        return Lit(sub_coord, ['xt'])

    def sort_limites(self):
        """
        Trie les limites par Xt croissants
        Étape nécessaire pour chercher les limites communes entre deux
          profils plus facilement
        """
        self.limites = self.limites.sort_index(by='Xt_profil')

    def calcul_nb_pts_inter(self, other, pas_long):
        """
        Calcul le pas sur la distance minimal entre les deux profils
        /!\ Cette valeur englobe les points de bord (début et fin)
        self et other sont commutables
        """
        dist_min_profil = self.geom.distance(other.geom)
        nb_pts_inter = ceil(dist_min_profil/pas_long) - 1
        return nb_pts_inter

    def export_profil_csv(self, outfile_profils, first, sep, digits):
        """
        Exporter les coordonnées dans un fichier CSV
        first <bool>:
          - True: effacement du fichier et ajout de l'entête avant écriture du profil
          - False: ajout du profil en fin de fichier
        """
        coord = self.coord.array

        # Ajout de la colonne id_profil
        # np_id_profil = np.empty(len(coord), dtype=np.int) #FIXME: modifier aussi avec le format
        np_id_profil = np.empty(len(coord), dtype=np.float)
        np_id_profil.fill(self.id)
        coord = append_fields(coord, 'id_profil', np_id_profil, usemask=False)

        fmt = sep.join(["%.{}f".format(digits)] * (len(coord.dtype.fields)-1) + ["%.1f"]) #FIXME: normalement format int pour id_profil
        if first:
            # Write header
            with open(outfile_profils, mode='w', newline='') as fileout:
                fileout.write(sep.join(coord.dtype.names) + '\n')
        with open(outfile_profils, mode='ab') as fileout:
            np.savetxt(fileout, coord, delimiter=sep, fmt=fmt)

    def export_limites_csv(self, outfile_limites, first, sep, digits):
        """
        Exporter les limites dans un fichier CSV
        """
        mode = 'w' if first else 'a'
        limites = deepcopy(self.limites)
        # Ajout des identifiants comme colonne
        limites['id_profil'] = self.id
        limites['id_ligne'] = limites.index
        limites.to_csv(outfile_limites, header=first, mode=mode, sep=sep, index=False, float_format="%.{}f".format(digits)) #FIXME: le format est plutôt pour id_profil aussi

    def change_coord(self, array):
        self.coord.array = array
        self.coord.compute_Xt()
        self.nb_points = len(self.coord.array)
        self.geom = self.coord.convert_as_linestring()


class SuiteProfilsTravers:
    """
    SuiteProfilsTravers: ensemble de profils en travers

    ### Attributs
    - suite <list <ProfilTravers>>

    ### Méthodes
    - find_and_add_limits
    - calculer_dist_proj_axe
    """
    def __init__(self, i3s_path, label, value_as_id=False):
        # super(SuiteProfilsTravers, self).__init__([])
        self.suite = []
        with BlueKenueRead_i3s(i3s_path) as in_i3s:
            in_i3s.read_header()
            for i, (value, linestring) in enumerate(in_i3s.iter_on_polylines()):
                id = value if value_as_id else i
                self.suite.append(ProfilTravers(id, list(linestring.coords), label))

    def __add__(self, other):
        newsuite = deepcopy(self)
        newsuite.suite = self.suite + other.suite
        return newsuite

    def __getitem__(self, index):
        if isinstance(index, slice):
            newsuite = deepcopy(self)
            newsuite.suite = self.suite[index]
            return newsuite
        else:
            return self.suite[index]

    def __len__(self):
        return len(self.suite)

    def __repr__(self):
        return [profil for profil in self.suite].__repr__()

    def find_and_add_limits(self, lignes_contraintes, dist_max):
        """
        lignes_contraintes <list <LigneContrainte>>
        """
        for i, profil_travers in enumerate(self):
            for ligne_contrainte in lignes_contraintes:
                profil_travers.find_and_add_limit(ligne_contrainte, dist_max)
            profil_travers.sort_limites()
            limits = list(profil_travers.limites.index)  # only to print

            print("> {}".format(profil_travers))
            print("{} limites trouvées avec les lignes {}".format(len(limits), limits))

            # Exports en CSV
            # profil_travers.export_profil_csv(args.outfile_profils, i, first_profil, sep, DIGITS)
            # profil_travers.export_limites_csv(args.outfile_limites, i, first_profil, sep, DIGITS)
            # if first_profil:
            #     first_profil = False

    def check_intersections(self):
        intersections = get_intersections([profil.geom for profil in self])
        if intersections:
            print("Les intersections suivantes sont trouvées")
            for (i, j) in intersections:
                print("- entre '{}' et '{}'".format(self[i], self[j]))
            sys.exit("ERREUR: Des profils s'intersectent")

    def calculer_dist_proj_axe(self, axe_geom):
        """
        Calculer la distance projetée sur l'axe
        axe_geom <LineString>
        /!\ Orientation de l'axe
        """
        for profil in self:
            profil_geom = profil.geom
            if profil_geom.intersects(axe_geom):
                # print(profil)
                # print(profil.coord.array)
                intersection = profil_geom.intersection(axe_geom)
                if isinstance(intersection, Point):
                    profil.dist_proj_axe = axe_geom.project(intersection)
                else:
                    sys.exit("L'intersection entre le '{}' et l'axe hydraulique n'est pas un point unique".format(profil))
            else:
                print(list(profil.geom.coords))
                sys.exit("ERREUR: Le '{}' n'intersection pas l'axe".format(profil))

    # def concat(self, other)


class LigneContrainte:
    """
    LigneContrainte: polyligne ouverte permettant de distinguer les lits

    ### Attributs
    - id <integer>: identifiant unique (numérotation automatique commençant à 0)
    - coord <2D-array float>: coordonnées ['X', 'Y'] des points de la ligne
    - geom <LineString>: objet géométrique

    ### Méthodes
    - coord_sampling_along_line
    """
    def __init__(self, id, coord):
        """
        Création d'un objet LigneContrainte à partir des coordoonées de la ligne
        coord: liste de tuples (de taille 2)
          Ex: [(x1,y1), (x2,y2), ..., (xn,yn)]
        """
        self.id = id
        self.nb_points = len(coord)
        self.coord = np.array(coord)
        self.geom = LineString(coord)

    def __repr__(self):
        return "Ligne de contrainte n°{} ({} points)".format(self.id, self.nb_points)

    @staticmethod
    def get_lines_from_i2s(i2s_path):
        """
        Extraction d'objects LineString à partir d'un fichier i2s
        Info: value is ignored
        """
        lines = []
        with BlueKenueRead_i2s(i2s_path) as in_i2s:
            in_i2s.read_header()
            for i, (value, linestring) in enumerate(in_i2s.iter_on_polylines()):
                lines.append(LigneContrainte(i, list(linestring.coords)))
        return lines

    def coord_sampling_along_line(self, Xp1, Xp2, Xp_adm_int):
        """
        Extraire les coordonnées des points interpolés le long de la ligne entre les deux abscisses
        Xp_adm_int <1D-array>: répartition des points entre Xp1 et Xp2 (valeurs de 0 à 1, avec un parcours de Xp1 à Xp2)
        Xp1, Xp2 <float>: abscisses de début et de fin
        """
        nb_pts_long = len(Xp_adm_int)

        # Construction de la liste des abscisses cibles (en mètres)
        Xp = (1 - Xp_adm_int)*Xp1 + Xp_adm_int*Xp2

        coord_int = []
        for Xp_cur in Xp:
            point = self.geom.interpolate(Xp_cur)
            coord_int.append(point.coords[0])
        np_coord_int = np.array(coord_int, dtype=float_vars(['X', 'Y']))

        return np_coord_int


class Lit(Coord):
    """
    Lit(s): morceau de profil borné par des limites
    ~> Hérite de la classe <Coord>

    ### Attributs
    - coord
    - nb_pts_trans

    ### Méthodes
    Toutes ces méthodes retournent un array avec les colonnes VARS2INT
    - interp_along_lit_auto
    - interp_along_lit
    - interp_inter_lineaire
    """
    VARS2INT = ['X', 'Y', 'Z']

    def interp_along_lit_auto(self, pas_trans):
        nb_pts_trans = ceil((self.array['Xt'][-1] - self.array['Xt'][0])/pas_trans) + 1
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_trans)
        return self.interp_along_lit(Xt_adm_list)

    def interp_along_lit(self, Xt_adm_list):
        """
        Echantionnage du lit (interpolation sur vecteur adimensionné de Xt)
        Xt_adm_list <1D-array float>: entre 0 et 1
        """
        nb_pts_trans = len(Xt_adm_list)
        array_ech = np.empty(nb_pts_trans, dtype=float_vars(Lit.VARS2INT))

        for label in Lit.VARS2INT:
            array_ech[label] = np.interp(Xt_adm_list, self.array['xt'], self.array[label])
        return array_ech

    def interp_inter_lineaire(self, other, coeff, nb_pts_trans):
        """
        Interpolation linéaire entre deux lits
        other <ProfilTravers>:
        coeff <float>: pondération entre self et other (0=self, 1=other)
        """
        # nb_pts_trans = ceil((1-coeff)*self.nb_pts_trans + coeff*other.nb_pts_trans)  #DEBUG: devrait dépendre de la "longueur" du profil aussi!
        Xt_adm_list = np.linspace(0., 1., num=nb_pts_trans)

        array_1 = self.interp_along_lit(Xt_adm_list)
        array_2 = other.interp_along_lit(Xt_adm_list)
        array = np.empty(nb_pts_trans, dtype=float_vars(Lit.VARS2INT))

        for label in Lit.VARS2INT:
            array[label] = (1-coeff)*array_1[label] + coeff*array_2[label]
        return array


class MeshConstructor:
    """
    MeshConstructor: données pour construire un maillage

    ### Attributs
    - points
    - i_pt (curseur pour répérer l'avancement)
    - segments

    ### Méthodes
    - add_points
    - add_segments
    - add_segments_from_node_list
    - export_as_dict

    /!\ Les éléments sont numérotés à partir de 0
      (dans le tableau segments)
    """
    def __init__(self):
        self.points = np.empty(0, dtype=float_vars(['X', 'Y', 'Z']))
        self.i_pt = -1
        self.segments = np.empty([0, 2], dtype=np.int)

    def add_points(self, coord):
        """
        Ajouter des sommets/noeuds
        coord <2D-array>: avec les colonnes ['X', 'Y', 'Z']
        """
        self.i_pt += len(coord)
        self.points = np.hstack((self.points, coord))

    def add_segments(self, seg):
        """
        Ajouter des segments à contraintre
        seg <2D-array int>: série de couples de sommets
        """
        self.segments = np.hstack((self.segments, seg))

    def add_segments_from_node_list(self, node_list):
        """
        Ajoute les sommets à partir d'une liste de noeuds successifs
        node_list <1D-array int>: série de noeuds
        """
        new_segments = np.column_stack((node_list[:-1], node_list[1:]))
        self.segments = np.vstack((self.segments, new_segments))

    def export_as_dict(self):
        """
        Exporter les données pour triangle.triangulate
        """
        return {'vertices': np.array(np.column_stack((self.points['X'], self.points['Y']))),
                'segments': self.segments}

    def debug(self, out_semis='debug_out_semis.xyz', out_i2s='debug_out_semis.i2s'):
        with open(out_semis, mode='wb') as fileout:
            np.savetxt(fileout, self.points, delimiter=' ', fmt='%.4f')

        with BlueKenueWrite_i2s(out_i2s, True) as i2s:
            i2s.write_header()
            for i, (n1, n2) in enumerate(self.segments):
                linestring = LineString([self.points[n1], self.points[n2]])
                i2s.write_polyline(linestring, i)
