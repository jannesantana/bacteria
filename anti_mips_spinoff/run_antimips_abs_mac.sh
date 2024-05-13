labels=()
values=()

while IFS=' ' read -r label value; do
    labels+=("$label")
    values+=("$value")
done < "$1"

# Execute with parameters
./antimips \
    "${values[0]}" \
    "${values[1]}" \
    "${values[2]}" \
    "${values[3]}" \
    "${values[4]}" \
    "${values[5]}" \
    "${values[6]}" \
    "${values[7]}"

exit 0