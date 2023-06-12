def parse_presets(preset) {
    // Protocol presets
    presets = [ "Zymo-Seq_miRNA" : [ "clip_R1":1, "three_prime_clip_R1":0, "adapter":"TGGAATTCTCGGGTGCCAAGG", "protocol":"illumina" ], 
                "illumina" : [ "clip_R1":0, "three_prime_clip_R1":0, "adapter":"TGGAATTCTCGGGTGCCAAGG", "protocol":"illumina" ],
                "nextflex" : [ "clip_R1":4, "three_prime_clip_R1":4, "adapter":"TGGAATTCTCGGGTGCCAAGG", "protocol":"nextflex" ],
                "qiaseq"   : [ "clip_R1":0, "three_prime_clip_R1":0, "adapter":"AACTGTAGGCACCATCAAT", "protocol":"qiaseq" ],
                "cats"     : [ "clip_R1":3, "three_prime_clip_R1":0, "adapter":"AAAAAAAA", "protocol":"cats" ]
    ]
    if (!(preset in presets.keySet())) {
        exit 1, "The provided protocol '${preset}' is not supported. The supported protocols are ${presets.keySet().join(', ')}"
    }
    clip_R1             = presets[preset]["clip_R1"]
    three_prime_clip_R1 = presets[preset]["three_prime_clip_R1"]
    adapter             = presets[preset]["adapter"]
    protocol            = presets[preset]["protocol"]
    return [clip_R1, three_prime_clip_R1, adapter, protocol]
}
