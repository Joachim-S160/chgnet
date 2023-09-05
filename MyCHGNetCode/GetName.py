def get_molecule_name(file: str, extra_tag: str = ""):
    """
    Args:   
        file: file to get the name of the molecule from
        extra_tag: extra tag to end at the end of the molecule name e.g.: _1000K
    returns: 
        molecule name
    """
    
    input_string = file

    # Get the index of the last occurrence of '/'
    last_slash_index = input_string.rfind('/')

    # Get the index of the first occurrence of '.'
    dot_index = input_string.find('.')

    # Extract the substring between the last slash and the first dot (exclusive)
    desired_substring = input_string[last_slash_index + 11:dot_index]

    # Find the third underscore
    third_underscore_index = desired_substring.find('_')

    if third_underscore_index != -1:
        desired_substring = desired_substring[:third_underscore_index]

    return desired_substring + extra_tag 

