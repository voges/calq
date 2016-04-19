int expand(str_t *exp, const char * const cigar, const char * const seq)
{
    size_t cigar_idx = 0;
    size_t cigar_len = strlen(cigar);
    size_t op_len = 0; // length of current CIGAR operation
    size_t seq_idx = 0;

    for (cigar_idx = 0; cigar_idx < cigar_len; cigar_idx++) {
        if (isdigit(cigar[cigar_idx])) {
            op_len = op_len * 10 + (size_t)cigar[cigar_idx] - (size_t)'0';
            continue;
        }

        size_t i = 0;
        switch (cigar[cigar_idx]) {
        case 'M':
        case '=':
        case 'X':
            // add matching part to expanded sequence
            str_append_cstrn(exp, &seq[seq_idx], op_len);
            seq_idx += op_len;
            break;
        case 'I':
        case 'S':
            seq_idx += op_len; // skip inserted part
            break;
        case 'D':
        case 'N':
            // inflate expanded sequence
            for (i = 0; i < op_len; i++) str_append_char(exp, 'D');
            break;
        case 'H':
        case 'P':
            break; // these have been clipped
        default: 
            cq_err("Bad CIGAR string: %s\n", cigar);
            return CQ_FAILURE;
        }

        op_len = 0;
    }

    return CQ_SUCCESS;
}