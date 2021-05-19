#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ValveMaskconverter.py
"""
create column ValveMask_label depending on value of ValveMask, Delta_18_16, and  Delta_D_H

DI  : δ18O(Delta_18_16):  -7.78 ± 0.01 permil; δD(Delta_D_H):  -50.38 ± 0.02 permil
GSM1: δ18O(Delta_18_16): -33.07 ± 0.02 permil; δD(Delta_D_H): -262.95 ± 0.04 permil
"""
# --- import -----------------------------------
# import from standard lib
import pandas as pd
from pathlib import Path
import yaml
from pprint import pformat
import argparse

# https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
# --- module's variable ------------------------
# global std


# ----------------------------------------------
class Standard(object):
    """ """
    def __init__(self, fpath_):
        """ initialise generic standard object
        :param fpath_: yaml file with parameters of standard
        """

        self._fpath = fpath_
        # list of attribute
        self._attr = ['delta', 'scale_factor', 'add_offset']

        #
        try:
            # read parameters configuration file yaml
            with open(fpath_, 'r') as stream:
                try:
                    self._param = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)
        except Exception:
            raise Exception(f"Something goes wrong when loading parameters file -{self._fpath}-.")

        # list all standard known
        self._set = {*self._param['default']['d18O']['std'],
                     *self._param['default']['dD']['std']
                    }
        self._isokey = self._param['default'].keys()
        #
        self.iso = {}
        # get default value for each
        self.get('default')

    def _unpacked(self):
        """
        'unpacked' standard gas concentration

        packed_value = (unpacked_value - add_offset) / scale_factor
        unpacked_value = (packed_value * scale_factor) + add_offset
        """

        try:
            for dx in self._isokey:
                scale_factor = self.iso[dx]['scale_factor']
                add_offset = self.iso[dx]['add_offset']
                for key in self._set:
                    _ = self.iso[dx][key]
                    _['unpacked'] = (_['packed'] * scale_factor) + add_offset
        except Exception:
            raise Exception(f"Something goes wrong unpacking data")


    def get(self, instr_):
        """ get parameters of instrument instr_

        overwrite default value if instrument value exist
        """
        try:
            _ = self._param[instr_]
            for dx in self._isokey:
                if dx in _:
                    # initialise dictionary if need be
                    if dx not in self.iso:
                        self.iso[dx] = {}
                    for key in self._attr:
                        if key in _[dx] and bool(_[dx][key]):
                            self.iso[dx][key] = _[dx].get(key)

                    if 'std' in _[dx] and bool(_[dx]['std']):
                        for key in self._set:
                            if key not in self.iso[dx]:
                                self.iso[dx][key] = {}
                            if bool(_[dx]['std'][key]):
                                # assume to be packed value
                                self.iso[dx][key]['packed'] = _[dx]['std'].get(key)

        except KeyError:
            raise KeyError(f"unknown instrument -{instr_}-. Check parameter file {self._fpath}")
        except Exception:
            raise Exception(f"Something goes wrong")

        self._unpacked()

    def range(self, name_, std_):
        """
        name_ ['d180', 'dD']
        std_ ['DI', 'GSM1']
        """
        try:
            d = self.iso[name_]
            try:
                unpacked_value = self.iso[name_][std_]['unpacked']

                mi = unpacked_value - self.iso[name_]['delta']
                ma = unpacked_value + self.iso[name_]['delta']
                return mi, ma
            except KeyError:
                raise KeyError(f"Uknown standard name -{name_}-")
            except Exception as e:
                raise Exception(f"{e}")
        except KeyError:
            raise KeyError(f"Uknown isotope name -{name_}-")
        except Exception as e:
            raise Exception(f"{e}")

    def label(self, d18, dD):
        """
        given concentration of 'Delta_18_16'(d18O) and  'Delta_D_H'(dD)
        return standard name
        """
        std_name = None
        for _ in self._set:
            d18_min, d18_max = self.range('d18O', _)
            dD_min, dD_max = self.range('dD', _)
            if d18_min <= d18 <= d18_max and dD_min <= dD <= dD_max:
                std_name = _
                break

        if std_name is None:
            std_name = "standard_unmatch"

        return std_name

    def label_ValveMask(self, row):
        """
        given value from ValveMask and Delta_18_16,  Delta_D_H,
        return label
        """
        if row['ValveMask'] == 0:
            return 'measurement'
        elif row['ValveMask'] == 6:
            return self.label(row['Delta_18_16_median'], row['Delta_D_H_median'])
        else:
            return 'unknown'


def add_ValveMask_label(fin_, fout_, fparam_):
    """
    """
    fpath = Path(fin_)

    std = Standard(fparam_)
    # get instrument name from filename
    instrument= str(fpath.name).split('-')[0]
    std.get(instrument)

    # Read data from file 'filename.csv'
    # (in the same directory that your python process is based)
    # Control delimiters, rows, column names with read_csv (see later)
    df = pd.read_csv(fpath, sep='\s+')
    # add 'rank' column
    # consecutive lines with the same 'ValveMask' value, get the same rank
    df['rank'] = ((df.ValveMask != df.ValveMask.shift()).cumsum())
    # compute median for each group of same rank
    df['Delta_18_16_median'] = df['Delta_18_16'].groupby(df['rank']).transform('median')
    df['Delta_D_H_median'] = df['Delta_D_H'].groupby(df['rank']).transform('median')
    # add column 'label', value depending on standard range
    # df['ValveMask_label'] = df.apply (lambda row: label_ValveMask(row, std), axis=1)
    df['ValveMask_label'] = df.apply (lambda row: std.label_ValveMask(row), axis=1)
    # remove temporary column, previously added
    df.drop(['rank', 'Delta_18_16_median', 'Delta_D_H_median'], axis=1, inplace=True)

    # write file
    df.to_csv(fout_, index=False)
    # with open(fout_, 'w') as f:
    #    df.to_string(f, index=False, col_space=25)
    # df.to_csv(, index=False)


def _parse():
    """set up parameter from command line arguments
    """
    # define parser
    parser = argparse.ArgumentParser(
        prog="Quince_labeliseValveMask",
        description="add name of the reference gas use to re-calibrate the sensor"
    )

    # positional arguments
    parser.add_argument("input",
                        type=str,
                        help="input directory or filename. If directory, all files with suffix '.dat' will be process")
    # optional arguments
    parser.add_argument("-o","--output",
                        type=str,
                        help="output directory or filename",
                        dest='output'
                        )
    parser.add_argument("-p","--parameter",
                        type=str,
                        help="yaml parameter configuration file, default from package directory",
                        default='parameters.yaml',
                        dest='param'
                        )

    # parse arguments
    args = parser.parse_args()

    return args


def main():
    """ """
    # read command line arguments
    args = _parse()

    # check arguments file
    pathin = Path(args.input)
    if pathin.is_file():
        dpathin = pathin.parent
        fpathin = pathin.name
    elif pathin.is_dir():
        dpathin = pathin
        fpathin = '*.dat'
    else:
        raise FileNotFoundError(f"can not find input file or directory -{pathin}-")

    if args.output is None:
        # use name of the input
        dpathout = dpathin
        fpathout = None
    else:
        pathout = Path(args.output)
        if pathout.is_dir():
            dpathout = pathout
            fpathout = None
        else:
            # Warning file to be created
            dpathout = pathout.parent
            fpathout = pathout.name

    fparam = Path(args.param).absolute()
    if not fparam.is_file():
        raise FileNotFoundError(f"can not find yaml parameters configuration file -{fparam}-")

    #
    for fin in dpathin.glob(fpathin):

        if fpathout is None:
            fout = dpathout / fin.with_suffix('.csv').name
        else:
            fout = dpathout / fpathout
            # add csv suffix if none
            if fout.suffix == '':
                fout =  fout.with_suffix('.csv')

        # check fin_ != fout_
        if fin == fout:
            raise FileExistsError(f"overwrite input file -{fin}-")
        # check fout_ already exist
        if fout.is_file():
            raise FileExistsError(f"file already exists -{fout}-")

        # print(fin, fout, fparam)
        add_ValveMask_label(fin, fout, fparam)


if __name__ == '__main__':
    """ """
    main()
