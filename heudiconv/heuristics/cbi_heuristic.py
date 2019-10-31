# Heuristics file in use at the NYU's Center for Brain Imaging
# (Our data come from a Siemens Prisma 3T scanner).
# Author: Pablo Velasco
# Date: 10/31/2019

import os

# TO DO: What do we do about pairs of scans "original" + "normalize" (when normalizing
#        for the coil sensitivity profile, but also keeping the original?)
#
#        Maybe I should use ".lower()" ???

def create_key(template, outtype=('nii.gz','dicom'), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return (template, outtype, annotation_classes)


def find_PE_direction_from_protocol_name( prot_name, default_dir_name='normal' ):
    PE_directions = ['AP','PA','RL','LR','rev']  # valid phase-encoding directions in the protocol name
    for direction in PE_directions:
        if ( ('_'+direction in prot_name) or
             ('-'+direction in prot_name) ):
            break    # we keep the current value of "direction"
        else:
            direction = default_dir_name    # fallback

    return direction



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
    t1 = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T1w')
    t2 = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T2w')
    pd_bias_body = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-biasBody_run-{item:02d}_PD')
    pd_bias_receive = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-biasReceive_run-{item:02d}_PD')
    # func:
    # For the functional runs, we want to have a different key for each task.  Since we don't know a-priori what the task
    #  names will be, we don't create the different keys here, but we'll create them dynamically, as needed.
    # diffusion:
    dwi = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_dwi')
    dwi_sbref = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_sbref')
    # fmaps:
    # For the SE-EPI field-map runs (for topup), both for fmri and dwi, we want to have a different key for
    #  each PE direction.  Rather than creating the different keys here, we'll create them dynamically, as needed.
    fmap_gre_mag = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-GRE_run-{item:02d}_magnitude')       # GRE fmap
    fmap_gre_phase = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-GRE_run-{item:02d}_phasediff')     #

    ###  DICOM only   ###
    # These are images we don't want for analysis, but we still want to keep a copy of
    #   the DICOMs in the 'sourcedata' folder.  We manage that by specifying the 'outtype'
    #   to be only 'dicom':
    # anat:
    t1_scout = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-Scout_run-{item:02d}_T1w', outtype = ('dicom',))
    t1_dicom = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T1w', outtype = ('dicom',))
    t2_dicom = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T2w', outtype = ('dicom',))
    # Misc:
    phoenix_doc = create_key('{bids_subject_session_dir}/misc/{bids_subject_session_prefix}_phoenix', outtype = ('dicom',))

    info = {t1:[], t2:[], pd_bias_body:[], pd_bias_receive:[],
            dwi:[], dwi_sbref:[],
            fmap_gre_mag:[], fmap_gre_phase:[],
            t1_scout:[], t1_dicom: [], t2_dicom: [], phoenix_doc:[]}

    for idx, s in enumerate(seqinfo):
        #pdb.set_trace()
        # s is a namedtuple with fields equal to the names of the columns
        # found in the dicominfo.tsv file

        acq = ''

        ###   T1w   ###
        # 1) Auto-Align Head scout (the original images, not the derived):
        #    look for "AA*Scout" or "AA*scout" in protocol_name:
        if ('scout' in s.protocol_name.lower()) and (s.is_derived == False):
            info[t1_scout].append({'item': s.series_id})
        # 2) High resolution T1w:
        # single volume, protocol name including T1, T1w, MPRAGE, MP-RAGE, MPR,...
        if ((s.dim4 == 1) and ( ('t1' in s.protocol_name.lower()) or
                                (('mp' in s.protocol_name.lower()) and ('rage' in s.protocol_name.lower())) or
                                ('mpr' in s.protocol_name.lower()) ) and
                              ('fl' in s.sequence_name)):
            # check the PE direction:
            direction = find_PE_direction_from_protocol_name( s.protocol_name, default_dir_name='' )
            acq = 'highres' + direction   # note: if direction is empty, aqc='highres'

            # If this image is NOT normalized, check if the previous or the following
            #   one has identical acquisition date and time.  If so, we'll keep only
            #   the normalized version, and for this one we keep just the DICOM.
            #   (Note: older versions of heudiconv don't include 'time'):
            # Otherwise, we extract it:
            if ( ('NORM' not in s.image_type) and
                 hasattr(s,'time') and
                 ( ( (idx+1 < len(seqinfo)) and (s.date == seqinfo[idx+1].date) and (s.time == seqinfo[idx+1].time) ) or
                     ( (idx > 0 ) and (s.date == seqinfo[idx-1].date) and (s.time == seqinfo[idx-1].time) ) ) ):
                info[t1_dicom].append({'item': s.series_id, 'acq': acq})
            else:
                info[t1].append({'item': s.series_id, 'acq': acq})

        # 3) FSE T1w:
        # single volume, series description includes TSE or FSE, protocol name includes T1, T1w
        if ((s.dim4 == 1) and ('t1' in s.protocol_name.lower()) and ('tse' in s.sequence_name)):
            # check the PE direction:
            direction = find_PE_direction_from_protocol_name( s.protocol_name, default_dir_name='' )
            acq = 'fse' + direction   # note: if direction is empty, aqc='fse'
            info[t1].append({'item': s.series_id, 'acq': acq})

        ###   T2w   ###
        # 1) Standard high-res T2w used for cortical segmentation:
        # single volume, protocol name including T2, T2w, TSE, SPACE, SPC:
        if (s.dim4 == 1) and (('T2' in s.protocol_name) or
                              ('tse' in s.protocol_name.lower()) or
                              ('space' in s.protocol_name.lower()) or
                              ('spc' in s.protocol_name.lower()) ):
            # check the PE direction:
            direction = find_PE_direction_from_protocol_name( s.protocol_name, default_dir_name='' )
            acq = 'highres' + direction   # note: if direction is empty, aqc='highres'

            # If this image is NOT normalized, check if the previous or the following
            #   one has identical acquisition date and time.  If so, we'll keep only
            #   the normalized version, and for this one we keep just the DICOM.
            #   (Note: older versions of heudiconv don't include 'time'):
            # Otherwise, we extract it:
            if ( ('NORM' not in s.image_type) and
                 hasattr(s,'time') and
                 ( ( (idx+1 < len(seqinfo)) and (s.date == seqinfo[idx+1].date) and (s.time == seqinfo[idx+1].time) ) or
                     ( (idx > 0 ) and (s.date == seqinfo[idx-1].date) and (s.time == seqinfo[idx-1].time) ) ) ):
                info[t2_dicom].append({'item': s.series_id, 'acq': acq})
            else:
                info[t2].append({'item': s.series_id, 'acq': acq})

        # 2) Fast Spin-Echo used as a localizer:
        # single volume, sequence name: 'h2d1' ('haste')"
        if (s.dim4 == 1) and ('h2d1' in s.sequence_name):
            info[t2_dicom].append({'item': s.series_id, 'acq': 'haste'})

        ###   PD   ###
        # BIAS images  (for coil sensitivity estimation) are typically PD-weighted
        if ('bias' in s.protocol_name.lower()) and ('tfl3d' in s.sequence_name):
            if ('body' in s.protocol_name.lower()):
                info[pd_bias_body].append({'item': s.series_id})
            else:
                info[pd_bias_receive].append({'item': s.series_id})

        ###   FUNCTIONAL   ###
        # We want to make sure the _SBRef and phase series (if present) are labeled
        #  the same as the main (magnitude) image.  So we only focus on the magnitude
        #  series (to exclude phase images) and more than 3 volumes (to exclude _SBRef)
        #  and then we search if the phase and/or _SBRef are present.
        if ((s.dim4 >= 4) and ('epfid2d' in s.sequence_name)
                          and (('M' in s.image_type) or ('FMRI' in s.image_type))
                          and (s.series_description[-6:].lower() != '_sbref')
                          and not ('DERIVED' in s.image_type)):

            ###   functional -- check PE direction   ###
            direction = find_PE_direction_from_protocol_name( s.protocol_name, default_dir_name='normal' )
            # we will write the 'direction' in under the 'acq-' tag.
            acq = direction

            ###   functional -- check task name   ###
            known_tasks=['rest', 'face', 'conflict', 'gamble', 'fix', 'films', 'inscapes']
            for task in known_tasks:
                if (task in s.protocol_name.lower()):
                    break    # we keep the current value for "task"
                else:
                    task = ''

            # if we don't find a known task, try finding the keyword "task":
            if ( task == '' ) and ('task' in s.protocol_name.lower()):
                # we want to capture what comes after "task", until the next
                #    dash ("-") or underscore ("_"):
                task = s.protocol_name.lower().split('task')[1]
                # remove any initial "-" or "_":
                if ((task[0] == '-') or (task[0] == '_')):
                    task = task[1:]
                # discard anything after the next "-" or "_":
                task = task.split('_')[0]
                task = task.split('-')[0]
                # remove any spaces we might have:
                task = task.replace(" ", "")
            else:
                task = 'TASK'    # fallback.  BIDS requires an alphanumeric string (no spaces...)

            # Below, for both magnitude and phase cases, we check if a template is defined
            #    with that specific name. That way we can assure the run number (run-01,
            #    run-02, ...) only apply for that specific task, not for just any
            #    task. E.g., if I have two tasks, with two different runs for each
            #    task name they should be called: task-foo_run-01, task-foo_run-02,
            #    task-boo_run-01,...    as opposed to task-boo_run-03 because there
            #    were 2 previous "task" runs.

            ###   functional -- is phase image present?   ###
            # At least for Siemens systems, if magnitude/phase was selected, the
            #  phase images come as a separate series immediatelly following the
            #  magnitude series.
            # (note: make sure you don't check beyond the number of elements in seqinfo...)

            if (idx+1 < len(seqinfo)) and ('P' in seqinfo[idx+1].image_type):
                # we have a magnitude/phase pair:

                # dictionary keys specific for this task type:
                mykey_mag = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_rec-magnitude_run-{item:02d}_bold" % task)
                if (mykey_mag in info):
                    info[mykey_mag].append({'item': s.series_id, 'acq': acq})
                else:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey_mag] = [{'item': s.series_id, 'acq': acq}]
                mykey_pha = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_rec-phase_run-{item:02d}_bold" % task)
                if (mykey_pha in info):
                    info[mykey_pha].append({'item': seqinfo[idx + 1].series_id, 'acq': acq})
                else:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey_pha] = [{'item': seqinfo[idx + 1].series_id, 'acq': acq}]

            else:
                # we only have a magnitude image

                # dictionary key specific for this task type:
                mykey = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_run-{item:02d}_bold" % task)
                if (mykey in info):
                    info[mykey].append({'item': s.series_id, 'acq': acq})
                else:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey] = [{ 'item': s.series_id, 'acq': acq }]

            ###   SB REF   ###
            # here, within the functional run code, check to see if the
            #  previous run protocol name ended in _SBREF, to assign the
            #  same task name and --if possible-- same run number.
            if (idx > 0) and ('_sbref' in seqinfo[idx - 1].series_description.lower()):
                mykey_sb = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_run-{item:02d}_sbref" % task)
                try:
                    # check if "info" already has this key by trying to append to it.
                    info[mykey_sb].append({'item': seqinfo[idx - 1].series_id, 'acq': acq})
                except KeyError:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey_sb] = [{'item': seqinfo[idx - 1].series_id, 'acq': acq}]


        ###   FIELD MAPS   ###

        # A) Spin-Echo distortion maps
        #    to be used with topup:
        #    (note: we only look for magnitude images, because we'll catch the phase images below)
        #pdb.set_trace()
        if ((s.dim4 <= 3) and ('epse2d' in s.sequence_name)
                          and ('M' in s.image_type)
                          and (  ('dist'  in s.protocol_name.lower())
                              or ('map'   in s.protocol_name.lower())
                              or ('field' in s.protocol_name.lower()) )):

            if (s.series_description[-4:] != '_SBRef'):    # sbref from MB diffusion have epse2d in
                                                           #  sequence_name, so don't include them
                ###   SE distortion: -- check PE direction   ###
                direction = find_PE_direction_from_protocol_name( s.protocol_name, default_dir_name='normal' )

                # Below, for both magnitude and phase cases, we check if a template is defined
                #    with that specific direction tag. That way we can assure the run number
                #    (run-01, run-02, ...) only applies for that specific direction tag, not
                #    for just any fmap. E.g., if I have two fmaps, with opposed directions
                #    they should be called: dir-AP_run-01 and dir-PA_run-01, as opposed to
                #    dir-AP_run-01 and dir-PA_run-02.  So each direction tag should be assigned
                #    to its own template.

                ###   SE-fmap -- is phase image present?   ###
                # At least for Siemens systems, if magnitude/phase was selected, the
                #  phase images come as a separate series immediatelly following the
                #  magnitude series.
                # (note: make sure you don't check beyond the number of elements in seqinfo...)

                if (idx+1 < len(seqinfo)) and ('P' in seqinfo[idx+1].image_type):
                    # we have a magnitude/phase pair:

                    # dictionary keys specific for this SE-fmap direction:
                    mykey_mag = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_rec-magnitude_dir-%s_run-{item:02d}_epi' % direction)
                    mykey_pha = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_rec-phase_dir-%s_run-{item:02d}_epi' % direction)
                    try:
                        # check if "info" already has this key by trying to append to it.
                        info[mykey_mag].append({'item': s.series_id})
                        info[mykey_pha].append({'item': seqinfo[idx + 1].series_id})
                    except KeyError:
                        # if it doesn't, add this key, specifying the first item:
                        info[mykey_mag] = [{'item': s.series_id}]
                        info[mykey_pha] = [{'item': seqinfo[idx + 1].series_id}]

                else:
                    # we only have a magnitude image

                    # dictionary key specific for this SE-fmap direction:
                    mykey = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_dir-%s_run-{item:02d}_epi' % direction)
                    try:
                        # check if "info" already has this key by trying to append to it.
                        info[mykey].append({'item': s.series_id})
                    except KeyError:
                        # if it doesn't, add this key, specifying the first item:
                        info[mykey] = [{ 'item': s.series_id}]



        # B) GRE fmap:
        if ('fm2d' in s.sequence_name) and (('fieldmap' in s.protocol_name.lower()) or
                                                ('field_map' in s.protocol_name.lower())):
#            pdb.set_trace()
            if (s.image_type[2] == 'M'):
                # magnitude image
                # For now, we don't need to separate according to TEs, because
                #   dcm2niix seems to create 2 NIfTI files appending "1" and "2"
                #   at the end, which is what BIDS specifies.
                # If this changes, we would have to do something, like s.TE
                info[fmap_gre_mag].append({'item': s.series_id})
            elif (s.image_type[2] == 'P'):
                # phase image:
                info[fmap_gre_phase].append({'item': s.series_id})


        ###   DIFFUSION   ###

        # We could also check: (s.image_type[2] == 'DIFFUSION')
        if ('ep_b' in s.sequence_name):       # Siemens product diffusion sequence:

            # This is not very rigorous, but in general, diffusion runs will
            #   have more than a couple of volumes.  I'll use 5, to be safe:
            if ( s.dim4 >= 5 ):
                # this is a standard diffusion acquisition
                acq = str(s.dim4)+'vols'
                info[dwi].append({'item': s.series_id, 'acq': acq})

                # check to see if the previous run is a SBREF:
                if ( (idx > 0) and
                         (seqinfo[idx - 1].series_description[-4:] == '_SBRef') and
                         ('epse2d' in seqinfo[idx - 1].sequence_name) ):
                    info[dwi_sbref].append({'item': s.series_id, 'acq': acq})

            else:
                # this is a fmap for diffusion.

                # TO-DO: for now, just ignore the _sbref image (if present, it would
                #        be the previous run (seqinfo[idx - 1].series_description[-4:] == '_SBRef')
                #        because BIDS doesn't allow them.

                # Get the PE direction, for topup:
                direction = find_PE_direction_from_protocol_name( s.protocol_name, default_dir_name='normal' )

                # dictionary key specific for this SE-fmap direction:
                mykey = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-%s_run-{item:02d}_epi' % direction)
                try:
                    # check if "info" already has this key by trying to append to it.
                    info[mykey].append({'item': s.series_id})
                except KeyError:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey] = [{ 'item': s.series_id}]

                if ( (idx > 0) and
                        (seqinfo[idx - 1].series_description[-4:] == '_SBRef') ):
                    mykey = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-%s_run-{item:02d}_sbref' % direction)
                    try:
                        # check if "info" already has this key by trying to append to it.
                        info[mykey].append({'item': seqinfo[idx - 1].series_id})
                    except KeyError:
                        # if it doesn't, add this key, specifying the first item:
                        info[mykey] = [{ 'item': seqinfo[idx - 1].series_id}]

        ###   PHOENIX FILE   ###

        if ('PhoenixZIPReport' in s.series_description) and (s.image_type[3] == 'CSA REPORT'):       #
            info[phoenix_doc].append({'item': s.series_id})


    #pdb.set_trace()
    return info
