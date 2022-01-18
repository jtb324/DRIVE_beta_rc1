import os


def listener(que_object, output: str, header: str):
    """
       Parameters
       __________
       que_object
           This parameter is the que object that is created earlier in
           the script. It is an object that the program can add data to
           in order to prevent data races

       output : str
           This parameter list the output path for the file

       header : str
           This parameter contains a string for the header row of the file

       """

    # opening the output file to write to
    with open(output, "a+") as output_file:

        # checking if the file size is zero
        if os.path.getsize(output) == 0:

            output_file.write(header)

        while 1:

            m = que_object.get()

            if m == "kill":

                break

            output_file.write(m)
            output_file.flush()
