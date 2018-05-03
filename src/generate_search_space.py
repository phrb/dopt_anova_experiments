import csv

yaml_path = ("../data/search_space.yaml")
yaml_file = open(yaml_path, "r")

# Drop first line
yaml_file.readline()

csv_line = []
headers = [ "elements_number",
            "y_component_number",
            "vector_length",
            "temporary_size",
            "vector_recompute",
            "load_overlap",
            "threads_number",
            "lws_y",
            "time_per_pixel"]

with open('../data/search_space.csv', 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(headers)

    for line in yaml_file:
        line = line.split()

        if line[0] == "?":
            csv_line = [line[2]]
        elif line[0] == ":":
            csv_line.append(line[1])
            writer.writerow(csv_line)
        else:
            csv_line.append(line[1])
