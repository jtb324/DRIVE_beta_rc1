import logging
from dataclasses import dataclass

@dataclass
class Logger():

    logging.basicConfig(
        level=logging.INFO,
        format='{asctime} {levelname:<8} {message}',
        style='{'
        )