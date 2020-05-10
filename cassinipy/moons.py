import csv
from datetime import datetime


class MoonFlyby:
    def __init__(self, flyby, instruments):
        self.flyby = flyby

        with open('data/moon_flyby_dates.csv', 'r', newline='') as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter=",")
            for row in csvreader:
                if row["flyby"] == flyby.upper():
                    self.date = datetime.strptime(row["date"], "%d-%b-%y")
                    self.moon = row["moon"]

    # def load_caps_data(self):

    # def closest_approach(self):


def moon_flybys(moon=""):
    """
    Returns a generator of moon flybys
    """
    with open('data/moon_flyby_dates.csv', 'r', newline='') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter=",")
        for row in csvreader:
            if moon == "":
                yield row.values()
            else:
                if row["moon"] == moon:
                    yield row.values()


def flyby_checker(flybydatetime):
    """
    Finds whether contains a moon flyby on a date, and returns moon/flyby or False
    """
    with open('data/moon_flyby_dates.csv', 'r', newline='') as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter=",")
        for row in csvreader:
            if row["date"] != "":
                flybydatetime = flybydatetime.strptime(row["date"], "%d-%b-%y")
                if data_datetime.date() == flybydatetime.date():
                    return row["moon"], row["flyby"]
    return False
