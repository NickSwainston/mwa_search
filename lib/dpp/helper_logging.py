import sys
import logging


class LoggerWriter(object):
    """
    I don't know how this works, I got it from stackoverflow
    It's used to redirect stderr/stdout to a logging instance
    """
    def __init__(self, writer):
        self._writer = writer
        self._msg = ''

    def write(self, message):
        self._msg = self._msg + message
        while '\n' in self._msg:
            pos = self._msg.find('\n')
            self._writer(self._msg[:pos])
            self._msg = self._msg[pos+1:]

    def flush(self):
        if self._msg != '':
            self._writer(self._msg)
            self._msg = ''


def initiate_logs(loglvl, outfile=None, writemode="a", stderr=True):
    """Initiates a logging instance"""
    logger = logging.getLogger()
    logger.setLevel(loglvl)
    formatter = logging.Formatter("%(asctime)s | %(filename)s | %(name)s | %(lineno)-4d | %(levelname)-9s || %(message)s",
                "%y/%m/%d %H:%M")
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    if outfile:
        fh = logging.FileHandler(outfile, writemode)
        fh.setLevel(loglvl)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    if stderr:
        sys.stderr = LoggerWriter(logger.error) # Redirect stderr to the log file