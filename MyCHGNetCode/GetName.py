def get_molecule_name(file: str):
    """
    Args: file
    returns: molecule name
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

    return desired_substring

