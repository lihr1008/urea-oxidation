class DictUtils:

    @classmethod
    def sort_count_dict(cls, count_dict, reverse: bool = False):
        sort_keys = sorted(count_dict.keys(), key=lambda x: count_dict[x], reverse=reverse)
        sort_dict = {}
        for k in sort_keys:
            sort_dict[k] = count_dict[k]
        return sort_dict


if __name__ == "__main__":
    pass
