import logging


def validateConfig(config: dict):
    """
    @description Validate the input config dict parsed from json

    @param config: dict
        config dict parsed from json
    """
    # NOTE:
    #   All config validation should be performed prior to calling the
    #   reduction function, following the "fail early" principle.

    # ------ #
    # Sample #
    # ------ #
    # --> multiple scattering correction
    # |
    # v Absorption correction
    #                   None  SampleOnly SampleAndContainer Mayers Carpenter
    # --------------------------------------------------------------------------
    # None               OK     ERR          ERR            ERR       ERR
    # SampleOnly         OK     OK           ERR            OK        WARN
    # SampleAndContainer OK     OK           OK             ERR       ERR
    # FullPaalmanPings   OK     OK           OK             ERR       ERR
    # Mayers             ERR    ERR          ERR            OK        ERR
    # Carpenter          ERR    ERR          ERR            ERR       WARN
    valid_absorption_methods = [
        "None",
        "SampleOnly",
        "SampleAndContainer",
        "FullPaalmanPings",
        "Mayers",
        "Carpenter",
    ]
    valid_multiple_scattering_methods = [
        "None",
        "SampleOnly",
        "SampleAndContainer",
        "Mayers",
        "Carpenter",
    ]
    sam_abs = config["Sample"].get("AbsorptionCorrection", None)
    sam_ms = config["Sample"].get("MultipleScatteringCorrection", None)
    if sam_abs is None:
        # missing AbsorptionCorrection section, which means no abs or ms
        # correction
        if sam_ms is None:
            logging.warning(
                "No AbsorptionCorrection or MultipleScatteringCorrection found."
            )
        else:
            raise ValueError(
                "AbsorptionCorrection section is missing but"
                "MultipleScatteringCorrection section is present"
            )
    else:
        sam_abs_type = sam_abs.get("Type", None)
        # quick check to make sure a known absorption correction method is used
        if str(sam_abs_type) not in valid_absorption_methods:
            raise ValueError(
                "AbsorptionCorrection.Type is invalid: {}".format(sam_abs_type)
            )
        # going through the logic map (see above)
        if sam_ms is None:
            logging.info(
                "The reduction will skip multiple scattering correction"
                "for sample."
            )
        else:
            sam_ms_type = sam_ms.get("Type", None)
            # quick check to make sure a known multiple scattering correction
            # method is used
            if str(sam_ms_type) not in valid_multiple_scattering_methods:
                raise ValueError(
                    "MultipleScatteringCorrection.Type is invalid: {}".format(
                        sam_ms_type
                    )
                )
            #
            if sam_abs_type is None:
                if sam_ms_type is not None:
                    raise ValueError(
                        "Any multiple scattering correction requires"
                        "corresponding absorption correction"
                    )
            elif sam_abs_type == "SampleOnly":
                if sam_ms_type == "SampleAndContainer":
                    raise ValueError(
                        f"ms_{sam_ms_type} is incompatible with"
                        f"abs_{sam_abs_type}."
                    )
                elif sam_ms_type == "Carpenter":
                    logging.warning(
                        "Carpenter multiple scattering correction is only"
                        "valid for Vanadium."
                    )
                else:
                    pass
            elif sam_abs_type in ["SampleAndContainer", "FullPaalmanPings"]:
                if str(sam_ms_type) not in ["None", "SampleAndContainer", "SampleOnly"]:
                    raise ValueError(
                        f"ms_{sam_ms_type} is incompatible with"
                        f"abs_{sam_abs_type}."
                    )
            elif sam_abs_type == "Mayers":
                if sam_ms_type != "Mayers":
                    raise ValueError(
                        f"ms_{sam_ms_type} is incompatible with"
                        f"abs_{sam_abs_type}."
                    )
            elif sam_abs_type == "Carpenter":
                if sam_ms_type != "Carpenter":
                    raise ValueError(
                        f"ms_{sam_ms_type} is incompatible with"
                        f"abs_{sam_abs_type}."
                    )
                else:
                    logging.warning(
                        "Carpenter multiple scattering correction is"
                        "only valid for Vanadium."
                    )
            else:
                logging.info(
                    f"[Sample] will use {sam_abs_type} for absorption,"
                    f"{sam_ms_type} for multiple scattering"
                )

    # ------------------------ #
    # Vanadium (normalization) #
    # ------------------------ #
    # --> multiple scattering correction
    # |
    # v Absorption correction
    #                    |   None	SampleOnly	 Mayers	Carpenter
    # ---------------------------------------------------------------
    # None               |    OK       ERR        ERR       ERR
    # SampleOnly         |    OK       OK         OK        OK
    # Mayers             |    OK       OK         OK        ERR
    # Carpenter          |    OK       OK         ERR       OK
    va_abs = config["Normalization"].get("AbsorptionCorrection", None)
    va_ms = config["Normalization"].get("MultipleScatteringCorrection", None)
    valid_absorption_methods = ["None", "SampleOnly", "Mayers", "Carpenter"]
    valid_multiple_scattering_methods = [
        "None", "SampleOnly", "Mayers", "Carpenter"]
    if va_abs is None:
        # missing AbsorptionCorrection section, which means no abs
        # or ms correction
        if va_ms is None:
            logging.warning(
                "Skip both absorption and multiple scattering correction."
            )
        else:
            raise ValueError(
                "Missing absorptionCorrection but found"
                "MultipleScatteringCorrection in Normalization"
            )
    else:
        va_abs_type = va_abs.get("Type", None)
        if str(va_abs_type) not in valid_absorption_methods:
            raise ValueError(
                "AbsorptionCorrection.Type is invalid: {}".format(va_abs_type)
            )
        #
        if va_ms is None:
            logging.info(
                "Will skip multiple scattering correction for normalization."
            )
        else:
            va_ms_type = va_ms.get("Type", None)
            if str(va_ms_type) not in valid_multiple_scattering_methods:
                raise ValueError(
                    "MultipleScatteringCorrection.Type is invalid: {}".format(
                        va_ms_type
                    )
                )
            #
            if va_abs_type is None:
                if va_ms_type is None:
                    logging.warning(
                        "skipping absorption and multiple scattering"
                        "for normalization."
                    )
                else:
                    raise ValueError(
                        "multiple scattering correction requires absorption"
                        "correction."
                    )
            elif va_abs_type == "Mayers":
                if va_ms_type == "Carpenter":
                    raise ValueError(
                        "Carpenter multiple scattering correction is"
                        "incompatible with Mayers absorption correction."
                    )
                else:
                    pass
            elif va_abs_type == "Carpenter":
                if va_ms_type == "Mayers":
                    raise ValueError(
                        "Mayers multiple scattering correction is"
                        "incompatible with Carpenter absorption correction."
                    )
                else:
                    pass
            else:
                logging.info(
                    f"[Normalization] will use {va_abs_type} for absorption"
                    f"and {va_ms_type} for multiple scattering"
                )
