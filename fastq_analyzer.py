def to_phred33(char: str) -> int:
    return ord(char) - 33


def sliding_window_triming(file_from: str, file_to: str, sliding_window_size=5, trimming_quality=30) -> None:
    reads_remaining = 0
    with open(file_from, 'r') as file, open(file_to, 'w') as file2:
        cur_line = 0
        num_trimmed_reads = 0
        file.seek(0)
        max_len_trimmed = 0
        min_len_trimmed = 500
        total_len_trimmed = 0
        buffer = []
        cnt = 0
        for line in file.readlines():
            line = line.strip()
            cur_line += 1
            if cur_line % 4 == 0:
                cnt += 1
                cutting_position = len(line)
                if cutting_position < sliding_window_size:
                    buffer = []
                    num_trimmed_reads += 1
                    continue
                sum_of_first = sum([to_phred33(line[i]) for i in range(sliding_window_size)])
                if sum_of_first <= trimming_quality * sliding_window_size:
                    buffer = []
                    num_trimmed_reads += 1
                    continue
                cur_phred = sum_of_first
                for i in range(len(line) - sliding_window_size):
                    cur_phred += to_phred33(line[i + sliding_window_size]) - to_phred33(line[i])
                    if cur_phred <= trimming_quality * sliding_window_size:
                        cutting_position = i + sliding_window_size
                        break
                while to_phred33(line[cutting_position - 1]) < trimming_quality and cutting_position > 1:
                    cutting_position -= 1
                if cutting_position < 1:
                    num_trimmed_reads += 1
                else:
                    total_len_trimmed += cutting_position
                    max_len_trimmed = max(max_len_trimmed, cutting_position)
                    min_len_trimmed = min(min_len_trimmed, cutting_position)
                    for i in range(len(buffer)):
                        good_string = buffer[i]
                        if i % 2 == 0:
                            len_pos = good_string.find("length")
                            good_string = f"{good_string[:len_pos]}length={cutting_position}"
                        elif i == 1:
                            good_string = good_string[:cutting_position]
                        else:
                            raise IOError("Буффер не корректно обрабатывается. Возможно, где-то не обновляется.")
                        file2.write(good_string + '\n')
                    file2.write(line[:cutting_position] + '\n')
                    reads_remaining += 1
                buffer = []
            else:
                buffer.append(line)
        print(f"Триммингу подверглось {num_trimmed_reads} прочтений")
        print(f"Максимальная длина прочтений в отфильтрованном файле: {max_len_trimmed}")
        print(f"Минимальная длина прочтений в отфильтрованном файле: {min_len_trimmed}")
        print(f"Средняя длина прочтений в отфильтрованном файле: {round(total_len_trimmed / reads_remaining)}")
        print(f"Оставшееся число прочтений в файле {file_to} после тримминга файла {file_from}: {reads_remaining}")
        print("------------------------------------")


file_name = input("Введите пожалуйста ТОЧНОЕ НАЗВАНИЕ файла для анализа:\n")
try:
    with open(file_name, 'r') as file:
        cur_line = 0
        max_len = 0
        min_len = 500
        sum_len = 0
        cnt = 0
        GC_content = 0
        phred33 = 0
        for line in file.readlines():
            cur_line += 1
            if cur_line % 4 == 1:
                line_splitted = line.split()
                for i in line_splitted:
                    if i.startswith("length"):
                        length = int(i[7:])
                        max_len = max(length, max_len)
                        min_len = min(length, min_len)
                        sum_len += length
                        break
                cnt += 1
            elif cur_line % 4 == 2:
                GC_content += line.count('G')
                GC_content += line.count('C')
            elif cur_line % 4 == 0:
                phred33 += to_phred33(line[9])

        print("Task 1")
        print(f"Общее число прочтений в файле равно {cnt}")
        print(f"Минимальная длина прочтения равна {min_len}")
        print(f"Средняя длина прочтения равна (округлил до целого числа) {round(sum_len / cnt)}")
        print(f"Максимальная длина прочтения равна {max_len}")

        print("------------------------------------")
        print("Task 2")
        print(f"GC-состав прочтений: {round(GC_content / sum_len * 100, 2)}%")

        print("------------------------------------")
        print("Task 3")
        print(f"Среднее значение качества по шкале Phred для всех прочтений для позиции 10:"
              f" {int(phred33 / cnt)}")

        print("------------------------------------")
        print("Task 4")
        file30 = "reads.fastq.window_5.qual30.txt"
        file60 = "reads.fastq.window_60.qual30.txt"
        sliding_window_triming("reads.fastq.txt", file30)
        sliding_window_triming(file30, file60, 60)
except FileNotFoundError:
    print("Некорректное название файла. Введите файл, который уже скачан или создан вами, и"
          " находящийся в той же директории, что и исходный код пролекта")