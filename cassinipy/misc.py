import calendar
from heliopy.util import config
import pathlib

def find_heliopy_download_dir():
    '''
    Find the install location of heliopy data directory
    '''
    with open(config.get_config_file(),"r") as file:
        for row in file:
            if row[:12] == "download_dir":
                heliopydatapath = pathlib.Path(row.rsplit("=")[1].rstrip())
    return heliopydatapath
find_heliopy_download_dir()

def toTimestamp(d):
    return calendar.timegm(d.timetuple())
