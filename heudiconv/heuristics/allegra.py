"""
Heuristics file in use at the NYU's Center for Brain Imaging
for old DWI data collected at the Allegra.
At the end of the classification, the function
`clean_up_unneeded_tags_from_info` simplifies the BIDS file
names

Author: Pablo Velasco
Date: 11/19/2020
"""

import os


DEFAULT_FIELDS = {
    # For the `dataset_description.json`:
    "Acknowledgements":
        "We thank Pablo Velasco and the rest of the NYU CBI (Center "
        "for Brain Imaging) personnel for preparing the BIDS "
        "dataset. "
        "TODO: whom else you want to acknowledge"
}


def create_key(subdir, file_suffix, outtype=('nii.gz', 'dicom'),
               annotation_classes=None, prefix=''):
    if not subdir:
        raise ValueError('subdir must be a valid format string')
    template = os.path.join(
        prefix,
        "{bids_subject_session_dir}",
        subdir,
        "{bids_subject_session_prefix}_%s" % file_suffix
    )
    return template, outtype, annotation_classes


def find_PE_direction_from_protocol_name(prot_name, default_dir_name='normal'):
    # valid phase-encoding directions in the protocol name
    PE_directions = ['AP','PA','RL','LR','rev']
    direction = default_dir_name
    for peDir in PE_directions:
        if (
            '_'+peDir in prot_name
            or '-'+peDir in prot_name
        ):
            direction = peDir
            break
    return direction


def add_series_to_info_dict(series_id, mykey, info, acq=''):
    """ adds a series to the 'info' dictionary """

    if info == None or mykey == '':
        return error

    if mykey in info:
        if acq == '':
            info[mykey].append({'item': series_id})
        else:
            info[mykey].append({'item': series_id, 'acq': acq})
    else:
        # if it isn't, add this key, specifying the first item:
        if acq == '':
            info[mykey] = [{'item': series_id}]
        else:
            info[mykey] = [{'item': series_id, 'acq': acq}]


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    
    allowed template fields - follow python string module: 
    
    item: index within category 
    subject: participant id 
    seqitem: run number during scanning
    subindex: sub index within group
    """
    
    ###   NIFTI and DICOM   ###
    # anat:
    t1 = create_key('anat','acq-{acq}_run-{item:02d}_T1w')
    t2 = create_key('anat','acq-{acq}_run-{item:02d}_T2w')

    # diffusion:
    dwi = create_key('dwi','acq-{acq}_run-{item:02d}_dwi')

    # fmaps:
    # For the SE-EPI field-map runs (for topup), both for fmri and dwi,
    # we want to have a different key for each PE direction.  Rather
    # than creating the different keys here, we'll create them
    # dynamically, as needed.

    ###  DICOM only   ###
    # These are images we don't want for analysis, but we still want to
    # keep a copy of the DICOMs in the 'sourcedata' folder.  We manage
    # that by specifying the 'outtype' to be only 'dicom':
    # anat:
    t1_scout = create_key('anat','acq-Scout_run-{item:02d}_T1w', outtype = ('dicom',))

    info = {t1:[], t2:[], dwi:[], t1_scout:[]}

    last_run = len(seqinfo)

    for s in seqinfo:
        # check the PE direction:
        direction = find_PE_direction_from_protocol_name(
            s.protocol_name, default_dir_name=''
        )
        # we don't need the original case of the protocol name:
        prot_name = s.protocol_name.lower()
        acq = ''

        if 'epse2d' in s.sequence_name:
            # Diffusion

            if '_shortte' in prot_name:
                # this is a fmap for diffusion.

                # Get the PE direction, for topup:
                run_dir = direction or 'normal'

                # dictionary key specific for this SE-fmap direction:
                mykey = create_key(
                    'fmap',
                    'acq-dwi_dir-%s_run-{item:02d}_epi' % run_dir
                )
                add_series_to_info_dict(s.series_id, mykey, info)

            else:
                # this is a standard diffusion acquisition
                nvols = prot_name.split('dwi_')[1].split('_b')[0]
                acq = str(nvols)+'vols'
                info[dwi].append({'item': s.series_id, 'acq': acq})


        else:
            pass

    return info
