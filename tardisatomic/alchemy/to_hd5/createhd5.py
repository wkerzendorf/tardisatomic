from abc import ABCMeta

import numpy as np

from tardis.atomic import plugins as atomic_plugins
from tardisatomic.alchemy import Ion, Atom, Level, Transition, TransitionType, \
    TransitionValue, TransitionValueType
from sqlalchemy.orm import aliased
import pandas as pd
import sys
import platform
import hashlib
import os
import gnupg
from time import gmtime, strftime

class BaseAtomicDatabase(object):
    __metaclass__ = ABCMeta

    atomic_db = None

    def _to_data_frame(self, sql_data):
        rec_data = [rec.__dict__ for rec in sql_data]
        return pd.DataFrame.from_records(rec_data)

    #ToDO:Finish that

    @atomic_db.setter
    def atomic_db(self, value):
        self.atomic_db =  value


class Atoms(atomic_plugins.Atoms, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])

    def load_sql(self):
        atomic_data = self.atomic_db.session.query(Atom).order_by(
            Atom.atomic_number).all()
        self.data = self._to_data_frame(atomic_data)


class Ions(atomic_plugins.Ions, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)


    def load_sql(self):
        ion_data = self.atomic_db.session.query(Ion, Atom).join("atom").values(
            Atom.atomic_number, Ion.ion_number,
            Ion.ionization_energy)
        self.data = self._to_data_frame(ion_data)


    def to_hdf(self, file_or_buf):
        raise NotImplementedError()


class Levels(atomic_plugins.Levels, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)


    def load_sql(self):
        level_data = self.atomic_db.session.query(Level, Ion, Atom).join(
            "ion", "atom").values(Atom.atomic_number, Ion.ion_number,
                                  Level.level_number, Level.g, Level.energy)
        self.dgnupgata = self._to_data_frame(level_data)


class Transitions(atomic_plugins.Transitions, BaseAtomicDatabase):
    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)
        self.transition_types = kwargs.get('transition_types', ['lines'])

    def load_sql(self):
        target_level = aliased(Level)
        source_level = aliased(Level)
        trans_a = aliased(Transition)
        trans_b = aliased(Transition)
        trans_data = \
            self.atomic_db.session.query(TransitionValue,
                                         TransitionValue.value,
                                         Transition,
                                         Transition.id.label(
                                             "transition_id"),
                                         target_level.id.label(
                                             "target_level_id"),
                                         source_level.id.label(
                                             "source_level_id"),
                                         target_level.level_number.label(
                                             "target_level_number"),
                                         source_level.level_number.label(
                                             "source_level_number"),
                                         source_level.g.label(
                                             'g_lower'),
                                         target_level.g.label(
                                             'g_upper'),
                                         TransitionType.name,
                                         TransitionValueType.name
                                         ).join(Transition).join(
                target_level,
                Transition.target_level_id == target_level.id).join(
                source_level,
                Transition.target_level_id == source_level.id
                ).join(
                TransitionType).filter(TransitionType.name == 'Line')

        raw_line_transitions = self._to_data_frame(trans_data)
        raw_line_transitions.drop('Transition', axis=1, inplace=True)
        raw_line_transitions.drop('TransitionValue', axis=1, inplace=True)
        raw_line_transitions_wl = raw_line_transitions.ix[raw_line_transitions[
                                                              'name'] ==
                                                          'wavelength']

        raw_line_transitions_wl.rename(columns={'value': 'wavelength'},
                                       inplace=True)
        raw_line_transitions_loggf = raw_line_transitions.ix[
            raw_line_transitions[
                'name'] == 'loggf']
        raw_line_transitions_loggf.rename(columns={'value': 'loggf'},
                                          inplace=True)
        self.line_transitions = pd.merge(raw_line_transitions_loggf[[
            'transition_id', 'loggf']],
                                         raw_line_transitions_wl,
                                         on='transition_id')

        self.line_transitions['f_ul'] = \
            self.line_transitions.apply(lambda x:
                                        np.power(10, x['loggf']) / x['g_upper'])
        self.line_transitions['f_lu'] = \
            self.line_transitions.apply(lambda x:
                                        np.power(10, x['loggf']) / x['g_lower'])
        self.data = self.line_transitions


class CreateHDF(object):
    def __init__(self, atom_db, hdf_file_or_buf ,exclude_species=[],
                 max_ionization_level=np.inf, \
                                                              transition_types=['lines']):
        """
        Creates the HDF file for TARDIS from the atomic database.

        :param atom_db: The AtomicDatabase object from
        tardisatomic.base containing all atomic data.

        :param exclude_species: Specifies the specie by atomic number which
        are exclude from the atomic data set. By default all species are
        included.

        :param max_ionization_level: Specifies the highest ionization level
        included in the atomic data set. All higher ionization levels are
        excluded. By default all ionization levels are included.

        :param transition_types: Specifies the transition types included in
        the atomic data set. Default: line
        """

        self.atomic_db = atom_db
        self.exclude_species = exclude_species
        self.max_ionization_level = max_ionization_level
        self.transition_types = transition_types

        atomic_plugins = self.__create_atomic_objects()



    def __create_atomic_objects(self):

        atomic_plugins ={x.hdf_name: x(exclude_species=self.exclude_species,
                                       max_ionization_level=
                                       self.max_ionization_level,
                                       transition_types =self.transition_types
                                       ) for x in
                         BaseAtomicDatabase.__subclasses__()}

        return atomic_plugins

    def __load_atomic_data(self,atomic_plugins):

        for key in atomic_plugins:
            try:
                atomic_plugins[key].load_sql()
            except:
                print("Unexpected error in {0}: {1} ".format(key,
                                                             sys.exc_info()[0]))



    def close_hdf(self,hdf_buf, anonymous=False):
        if not anonymous:
            uname = platform.uname()
            user = os.getusername()
        else:
            uname = 'anonymous'
            user = 'anonymous'
        ctime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        md5_hash = hashlib.md5()

        gpg = gnupg.GPG()
        gpg.encoding = 'utf-8'

        hdf_buf.put('metadata/system', uname)
        hdf_buf.put('metadata/username', user)
        hdf_buf.put('metadata/creation_time', ctime)

        for dataset in hdf_buf.values():
            md5_hash.update(dataset.value.data)

        data_to_sign = "TARDISATOMIC\n" "creation_time: " + ctime+ "\n" + \
            "System:\n" + \
                       uname +"\n"+ "User: " \
                                                                      "" + \
                       user + "\n" + "md5: " + md5_hash + "\n"
        signed_data = gpg.sign(data_to_sign)

        hdf_buf.put('metadata/gpg',  signed_data.data)
        hdf_buf.put('metadata/md5',md5_hash.hexdigest() )
        hdf_buf.close()













