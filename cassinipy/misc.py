import calendar
import pathlib

from heliopy.util import config


def find_heliopy_download_dir():
    """
    Find the install location of heliopy data directory
    """
    with open(config.get_config_file(), "r") as file:
        for row in file:
            if row[:12] == "download_dir":
                heliopydatapath = pathlib.Path(row.rsplit("=")[1].rstrip())
    return heliopydatapath


def toTimestamp(d):
    return calendar.timegm(d.timetuple())
