# This file containes the function that is used to determine if an individual contains a desired variant
def individualVariantID(row, variantList):
    '''this is a function for the inner loop that will search through each position in the row and when it encouters a one or a two it will add that to the idlist and then return so that the outer loop in the main script moves on to the next row.'''

    for i in range(6, (len(row)-1)):  # This for loop loops through the length of the recode row

        if row[i] == '1' or row[i] == '2':  # Checks to see if the element at row[i] is either 1 or 2

            # If it is 1 or 2 then it appends that id to totalVariantList
            variantList.append(row[1])

            return  # This return will then break out of the loop and move on to the next row
