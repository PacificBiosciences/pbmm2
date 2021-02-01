function sorttags() {
    for i in $("$SAMTOOLS" view $1); do echo -n "$i" | cut -f 1-11 | awk '{printf($0)}'; echo -n "$i" | cut -f 12- - | sort | tr '\n' '\t' | awk '{print "\t"$0}' ; done
}
