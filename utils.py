def combine_list_of_dictionaries(dlist):
    result_dict = {}
    for dic in dlist:
        for key in dic.keys():
            if not key in result_dict.keys():
                result_dict[key] = []
            result_dict[key].append(dic[key])

    return result_dict