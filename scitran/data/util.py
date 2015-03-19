# @author:  Gunnar Schaefer

import calendar
import datetime


def datetime_encoder(o):
    if isinstance(o, datetime.datetime):
        if o.utcoffset() is not None:
            o = o - o.utcoffset()
        return {"$date": int(calendar.timegm(o.timetuple()) * 1000 + o.microsecond / 1000)}
    raise TypeError(repr(o) + " is not JSON serializable")


def datetime_decoder(dct):
    if "$date" in dct:
        return datetime.datetime.utcfromtimestamp(float(dct["$date"]) / 1000.0)
    return dct
