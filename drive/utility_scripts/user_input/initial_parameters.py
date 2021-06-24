import typer

# This script will ask for the user to input certain initial parameters

def check_provided_integer(value: int) -> int:
    """callback function that will be used to make 
    sure the MIN_CM and the THREADS values falls 
    within a required range 
    Parameters
    __________
    value : int
        int that has to be greater than 0
    
    Returns
    _______
    int
        returns the integer value if it is valid. Otherwise it raises an error if the value was 
        not an integer or if the value is not greater than or equal to 0
    """
    if type(value) != int:
        raise typer.BadParameter(f"You passed a value of {value} as an argument. The value is expected to be an integer")
    if value < 0:
        raise typer.BadParameter(f"You passed a value of {value}, the program expects a positive value")
    return value

def check_provided_maf(value: float) -> float:
    """callback function that will be used to make sure the maf threshold falls within a required range 
    Parameters
    __________
    value : float
        float that falls between 0-0.5. This value 
        will default to 0.5 if the user doesn't 
        provide a value
    
    Returns
    _______
    float
        returns the float value if it is valid. Otherwise it raises an error
    """
    if not 0 <= value < 0.5:
        raise typer.BadParameter(f"The value {value} was passed as an argument. The program expected this parameters to be within the range [0,0.5)")
    return value


