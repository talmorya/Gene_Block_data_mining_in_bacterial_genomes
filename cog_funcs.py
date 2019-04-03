def get_cogs_info(cog_list):
    cog_dic = {}

    with open('/Users/Shani/PycharmProjects/Gene_Block_data_mining_in_bacterial_genomes/COG_INFO_TABLE.txt') as f:
        line = f.readline()

        # while there is a next line parse it
        while line:

            cog_name = line[3:7]
            if cog_name in cog_list:
                line_elements = line.split(';')
                cog_dic[cog_name] = line_elements[2:]

            line = f.readline()

    return cog_dic
