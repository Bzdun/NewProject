﻿CREATE TABLE [dbo].[AuthoringTable] (
    [ID]             UNIQUEIDENTIFIER CONSTRAINT [DF_AuthoringTable_ID] DEFAULT (newid()) NOT NULL,
    [AuthoringTable] NVARCHAR (100)   NOT NULL,
    [CreatedBy]      UNIQUEIDENTIFIER NULL,
    [ModifiedBy]     UNIQUEIDENTIFIER NULL,
    [CreatedOn]      SMALLDATETIME    CONSTRAINT [DF_AuthoringTable_CreatedOn] DEFAULT (CONVERT([smalldatetime],getdate(),(0))) NULL,
    [ModifiedOn]     SMALLDATETIME    CONSTRAINT [DF_AuthoringTable_ModifiedOn] DEFAULT (CONVERT([smalldatetime],getdate(),(0))) NULL,
    [IsDeleted]      BIT              NULL,
    [CurrVersion]    INT              NULL,
    CONSTRAINT [PK_AuthoringTable] PRIMARY KEY CLUSTERED ([ID] ASC),
    FOREIGN KEY ([CurrVersion]) REFERENCES [dbo].[CZVersion] ([VersionNumber]) ON DELETE NO ACTION ON UPDATE NO ACTION,
    CONSTRAINT [FK_AuthoringTable_User] FOREIGN KEY ([CreatedBy]) REFERENCES [dbo].[User] ([ID]) ON DELETE NO ACTION ON UPDATE NO ACTION,
    CONSTRAINT [FK_AuthoringTable_User1] FOREIGN KEY ([ModifiedBy]) REFERENCES [dbo].[User] ([ID]) ON DELETE NO ACTION ON UPDATE NO ACTION
);

