﻿








CREATE PROCEDURE [dbo].[GenerateBibliographyViewVersion]
@VersionNumber int
AS

INSERT INTO dbo.Staging_BibliographyView
SELECT *
FROM (
SELECT 
	Version_ExhibitReference.ExhibitID,
	Version_ExhibitReference.IsDeleted,
	Version_Reference.ID,
	Version_Reference.Title,
    Version_Reference.Authors,
    Version_Reference.BookChapters,
    Version_CitationType.CitationType,
    Version_Reference.PageNumbers,
    Version_Reference.Publication,
    Version_Reference.PublicationDates,
    Version_Reference.Source,
    Version_ExhibitReference.CurrVersion,
    CASE 
		WHEN (Version_ExhibitReference.IsDeleted = 1) OR 
			(Version_Reference.IsDeleted = 1) OR 
			(Version_Reference.IsVisible = 0) OR 
			(Version_CitationType.IsDeleted = 1) THEN 0
			ELSE 1 END AS [IncludeInVersion],
		NEWID() as NEW_ID
FROM Version_ExhibitReference
INNER JOIN Version_Reference ON Version_Reference.ID = Version_ExhibitReference.ReferenceID
			AND	Version_Reference.CurrVersion = (	SELECT MAX(CurrVersion) 
											FROM Version_Reference 
											WHERE Version_Reference.CurrVersion <= @VersionNumber
											AND Version_Reference.ID = Version_ExhibitReference.ReferenceID)
INNER JOIN Version_CitationType ON Version_CitationType.ID = Version_Reference.CitationTypeID
			AND	Version_CitationType.CurrVersion = (	SELECT MAX(CurrVersion) 
											FROM Version_CitationType 
											WHERE Version_CitationType.CurrVersion <= @VersionNumber
											AND Version_CitationType.ID = Version_Reference.CitationTypeID)
WHERE Version_ExhibitReference.CurrVersion = (	SELECT MAX(CurrVersion) 
								FROM Version_ExhibitReference as CurrView
								WHERE CurrView.CurrVersion <= @VersionNumber
								AND CurrView.ExhibitID = Version_ExhibitReference.ExhibitID)
) as tab
WHERE [IncludeInVersion] = 1







