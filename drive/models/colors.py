from enum import Enum, unique

@unique
class Colors(Enum):
    """A class to have multiple colors for warnings and errors and success."""
    RED: str = '\033[0;31m'
    YELLOW: str = '\033[1;33m'
    GREEN: str = '\033[0;32m'
    NOCOLOR: str = '\033[0m'    