
def has_ambiguous_bases(sequence):
    # Check if sequence contains lower-case letters, which usually indicate soft-masked bases
    ambiguous_bases = ['N', 'X', 'n', 'x']
    return any(base in ambiguous_bases for base in sequence)

def select_from_tsv(input_lines, selected_columns=None, schema=None, result_structure=tuple, default_column_type=str):
    """
    Parse selected_columns from a table represented in tab-separated-values format (TSV).

    For example, given an input_stream of lines

        "name\tage\theight\n",
        "tommy\t9\t120\n",
        "masha\t7\t122\n",
        ...

    the generator

        select_from_tsv(input_stream, selected_columns=["name", "height"])

    would yield

       ("tommy", "120"),
       ("masha", "122"),
       ...

    Note how in the above example all values returned are strings.

    Column type information may be specified as follows

        select_from_tsv(input_stream, selected_columns={"name": str, "height": float})

    yielding

        ("tommy", 120.0),
        ("masha", 122.0),
        ...

    If the first line in the input does not list all column headers, the schema must be provided as an additional argument

        select_from_tsv(input_stream, selected_columns=["name", "height"], schema={"name": str, "age": int, "height": float})

    Specifyng a schema argument means the first line of input should contain values, not columnd headers.

    Type information may be specified in either the selected_columns or the schema argument.

    Requesting result_structure=dict will change the output to

       {"name": "tommy", "height": 120.0},
       {"name": "masha", "height": 122.0},
       ...

    When a schema argument is provided, leaving selected_columns unspecified is equivalent to selecting all columns from the schema.

    At least one of the schema or selected_columns arguments must be specified, even if the schema could be inferred from the input headers.  This makes applications more robust.
    """
    assert schema != None or selected_columns != None, "at least one of these arguments must be specified"
    lines = strip_eol(input_lines)
    j = 0
    if schema == None:
        # the first line is expected to list all column names (headers)
        headers = next(lines).split('\t')  # pylint: disable=stop-iteration-return
        schema = headers
        j = 1
    if not isinstance(schema, dict):
        schema = {c: default_column_type for c in schema}
    headers = list(schema.keys())
    # Ensure "selected_columns" is an ordered dict of {column_name: column_type} and "schema" is a superdict of "columns".
    if selected_columns == None:
        # Return all columns
        selected_columns = dict(schema)
    if not isinstance(selected_columns, dict):
        selected_columns = {c: default_column_type for c in selected_columns}
    # Merge column type information from schema and selected_columns
    for c in schema:
        if schema[c] == default_column_type:
            schema[c] = selected_columns.get(c, schema[c])
    for c in selected_columns:
        if selected_columns[c] == default_column_type:
            selected_columns[c] = schema.get(c, selected_columns[c])
    column_indexes = []
    for c in selected_columns:
        assert c in headers, f"Column not found in table headers: {c}"
        column_indexes.append(headers.index(c))
        assert schema[c] == selected_columns[c], f"Conflicting types for column {c}."
    column_types = list(selected_columns.values())
    column_names = list(selected_columns.keys())
    for i, l in enumerate(lines):
        values = l.split('\t')
        if len(headers) != len(values):
            assert False, f"Line {i + j} has {len(values)} columns;  was expecting {len(headers)}."
        # type-convert and reorder values as specified in schema
        ordered_values = (ctype(values[ci]) for ci, ctype in zip(column_indexes, column_types))
        if result_structure in (tuple, list):
            yield result_structure(ordered_values)
        else: # dict
            yield result_structure((c, val) for c, val in zip(column_names, ordered_values))

def strip_eol(lines_iterable):
    for line in lines_iterable:
        yield line.rstrip('\n')